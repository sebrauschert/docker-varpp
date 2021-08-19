#!/usr/bin/env Rscript

# R script to rank variants in a patient .vcf file
# 

# Load the necessary libraries
library(optparse)

#=============================================================================================================================

# OPTPARSE options
option_list = list(
  make_option(c("-v", "--patient_vcf"),
              type="character",
              default=NULL,
              help="A patient derived .vcf file.",
              metavar="character"),
  make_option(c("-d", "--data"),
              type="character",
              default=NULL,
              help="Model data to predict with. Either GTEx or HCL.",
              metavar="character"),
  make_option(c("-p", "--hpo"),
              type="character",
              default=NULL,
              help="HPO term ID to predict with. This needs to be in the list of HPO terms that were pre-modelled.",
              metavar="character"),
  make_option(c("-c", "--custom_model"),
              type="character",
              default=NULL,
              help="Path to a pre-trained VARPP-RuleFit model in the form of a .rds object If this option is chosen,-p/--hpo option will be ignored.",
              metavar="character"),
  make_option(c("-o", "--output_path"),
              type="character",
              default="results/",
              help="Specify the path to the results. Defaults to the 'results' in the current directory.",
              metavar="character")
);

opt_parser <- OptionParser(option_list = option_list,
                           description = "\n varpp prediction script.",
                           epilogue = "Example:\n\n  ./predict.R -v <patient.vcf> -d <data: GTEx or HCL> -p <HPO term ID> -o <output/path> \n\n");

opt <- parse_args(opt_parser)

if (is.null(opt$patient_vcf)){
  print_help(opt_parser)
  stop("Missing patient .vcf file!\n", call.=FALSE)
}

patient_file_path  <- paste0('/mnt/', opt$patient_vcf)
output_file_path   <- paste0('/mnt/', opt$output_path, '/')
model_data         <- opt$data
hpo_term           <- opt$hpo
custom_model        <- opt$custom_model
#=============================================================================================================================
# Below are the fixed paths and the code
library(tidyverse)
library(vcfR)
library(varppRule)

# Set the path for the patient file and results file
genome_file_path  <- paste0("/mnt/genome-data/",opt$data) # this should be fixed either gtex or hcl (currently only gtex)
scratch           <- "/mnt/scratch/" # directoy where we save the intermediate files 
header_file       <- "/root/dependencies/" # This is part of the docker image

#================================================================================================================================
# Functions
#================================================================================================================================

range01 <- function(x){
  
  (x-min(x))/(max(x)-min(x))
  
}

# Prediction function based on the rules of the hpo terms.
# This fuction will return a list of variants that is potentially disease causing
# Might include CADD score for this.

# Logit to prop
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}


# the predict function for VARPP-RuleFit
predict.varppRule <- function(patient_data, 
                              results, 
                              predict=c("probability", "class"), 
                              type = c("GTEx","HCL")){
  
  attr(results, "class") <- "varppRule" # This is legacy, as all results have been created with an old attribute name
  
  # Prepare coefficients
  rules <- results$RuleFit$varimp[grep("rule",results$RuleFit$varimp$Variable),] %>%
    select(Description) %>%
    unlist() %>%
    as.character()
  
  coefficients <- results$RuleFit$varimp[grep("rule",results$RuleFit$varimp$Variable),] %>%
    select(Coefficient) %>%
    unlist() %>%
    as.numeric()
  
  
  
  # Create Rule variables based on the patient data
  for( i in 1:length(rules)){
    
    if(grepl('CADD_raw_rankscore >', rules[i])) {
      patient_data <-  cbind(patient_data, patient_data %>%
                               mutate(!!rules[i] := ifelse(eval(parse(text=rules[i])), 1, 0 )) %>%
                               select(!!rules[i]) )
      
    }else{
      
      patient_data <- cbind(patient_data,patient_data %>%
                              mutate(!!rules[i] := ifelse(eval(parse(text=rules[i])), 0, 1 )) %>%
                              select(!!rules[i]))
    }
    
  }
  
  # Prediction step
  features <- as.character(results$RuleFit$varimp$Variable)
  
  # Make sure the rule variables are called the same in the feature set as they were called above in the matrix
  features[grep("rule", results$RuleFit$varimp$Variable)] <-
    as.character(results$RuleFit$varimp$Description[grep("rule", results$RuleFit$varimp$Variable)])
  
  feature_vector <- results$RuleFit$varimp$Coefficient
  names(feature_vector) <- features
  
  # Here we return the probabilities for every variant to be pathogenic
  patient_data$Prediction <- logit2prob(rowSums(as.matrix(patient_data[,features]) %*%diag(feature_vector)))
  
  if(predict=="probability"){
    
    patient_data[,c("Chr", "Start", "End","Gene", "CADD_raw_rankscore", "Prediction")]
    
  }else{
    
    patient_data %>%
      mutate(Prediction = ifelse(Prediction > 0.5, "Pathogenic", "Benign")) %>%
      select(Chr, Start, End, Gene, CADD_raw_rankscore, Prediction)
  }
  
}

#================================================================================================================================

# Read the patient vcf file
patient_variants <- read.vcfR(patient_file_path)

# Bring the file into the right format (bedfile) and make sure the 'chr' is appended to the chromosome (vcf files not including this?)
as_tibble(getFIX(patient_variants)) %>%
  mutate(CHROM = paste0("chr", CHROM)) %>%
  mutate(Start = as.numeric(POS) - 1) %>%
  mutate(End = as.numeric(POS)) %>%
  select(CHROM, Start,End) -> dat

# Save the vcf so we can use bedtools
write_delim(dat, paste0(scratch,"patient.bed"), delim="\t", col_names=F)

#================================================================================================================================
## # Run bedtools sort and merge the patient variants with the data necessary for the prediction: GTEx or HCL
system(paste0("/root/miniconda3/bin/bedtools sort -i ", scratch, "patient.bed > ", scratch, "patient_sorted.bed;
cd ", genome_file_path," ;
ls *.bed.gz | parallel '/root/miniconda3/bin/bedtools  intersect -a ", scratch, "patient_sorted.bed -b {} -wb -sorted > ", scratch, "{.}_complete'"), wait=TRUE)

# Concatenate all the individual files to one result file
system(paste0("cat ", header_file, "header.txt > ", scratch,"mergelist.txt;
for i in $(seq 1 22) X Y; do cat ", scratch,"SORTED_FINAL_chr$i.bed_complete | cut -f1,2,3,8,9 --complement >>", output_file_path, "mergelist.txt; done
"), wait=TRUE)

# Remove the unnecessary files
#system(paste0('rm ',scratch, '*.bed_complete ', scratch, '*.bed'))

## # Read in the results file
predict <- read_delim(paste0(output_file_path, "mergelist.txt"), delim = "\t")
predict <- predict %>%
  mutate(CADD_raw_rankscore = CADD_2)


# Clean up all the unnecessary files and only retain the patient vcf for prediction
system(paste0("rm ", scratch, "*.bed_complete ", scratch, "patient.bed ", scratch, "patient_sorted.bed ", output_file_path, "mergelist.txt" ))

#================================================================================================================================
# Predict

if(is.null(custom_model)){
  
  model       <- readRDS(paste0('/mnt/models/',model_data,'/',hpo_term,'_RuleFit.rds')) # in case VARPP is run in the same session, the model can directly be assigned.

} else{
  
  model <- readRDS(custom_model)  
  
}
PREDICTIONS <- predict.varppRule(predict,
                                 model,  
                                 predict="probability")

#================================================================================================================================
# Print out header of results
cat("Predictions are:\n")
head(PREDICTIONS)

# Write the results
write_csv(PREDICTIONS, paste0(output_file_path,"predictions.csv"))
