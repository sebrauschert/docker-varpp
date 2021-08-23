#!/usr/bin/env Rscript

# Script to run VARPP-RuleFit and return a prediction
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
  make_option(c("-g", "--gene_list"),
              type="character",
              default=NULL,
              help="A comma seperated list of potentially pathogenic genes associated with the patient.(e.g. ABC,DEF,GHI)",
              metavar="character"),
  make_option(c("-t", "--type"),
              type="character",
              default=NULL,
              help="Run VARPP-RuleFit with either 'hcl', 'gtex' or 'custom' data",
              metavar="character"),
  make_option(c("-p", "--user_pathogenic"),
              type="character",
              default=NULL,
              help="A file with custom pathogenic data in the same format as the HCL or GTEx data.",
              metavar="character"),
  make_option(c("-b", "--user_benign"),
              type="character",
              default=NULL,
              help="A file with custom benign data in the same format as the HCL or GTEx data.",
              metavar="character"),
  make_option(c("-o", "--output_path"),
              type="character",
              default="results/",
              help="Specify the path to the results. Defaults to the 'results' in the current directory.",
              metavar="character"),
  make_option(c("-c", "--cores"),
              type="character",
              default=4,
              help="Maximum number of cores available. This will be distributed across the model",
              metavar="character"),
  make_option(c("-n", "--ntree"),
              type="character",
              default=2000,
              help="Number of trees in the model, defaults to 2000.",
              metavar="character"),
  make_option(c("-l", "--lasso"),
              type="character",
              default=100,
              help="Number of LASSO bootstrap rounds. Defaults to 100.",
              metavar="character"),
  make_option(c("-m", "--max_depth"),
              type="character",
              default=3,
              help="Maximum tree depth; defaults to 3..",
              metavar="character")
);

opt_parser <- OptionParser(option_list = option_list,
                           description = "\n VARPP-RuleFit.",
                           epilogue = "Example:\n\n  ./varppRuleAndPredict.R -v <path to patient vcf file> -g <gene list> -t <type: GTEx, HCL, custom> -p <User provided pathogenic data> -b <User provided benign data> -o <output/path> -c <Maximum number of cores available> -n <Number of trees, defaults to 2000> -l <Number of LASSO bootstrap rounds> -m <Tree depth>\n\n");

opt <- parse_args(opt_parser)

library(varppRule)
library(tidyverse)
library(vcfR)

genes       <- str_trim(unlist(str_split(opt$gene_list, ",")))
type        <- opt$type
user_patho  <- opt$user_pathogenic
user_benign <- opt$user_benign
output_path <- opt$output_path
cores       <- as.numeric(opt$cores)
ntree       <- as.numeric(opt$ntree)
patient_vcf <- opt$patient_vcf
scratch     <- '/mnt/scratch/'
max.depth   <- as.numeric(opt$max_depth)
lasso       <- as.numeric(opt$lasso)

# Run the model
model <- rule_fit(HPO_genes = genes,
                    HPO_term_name = "custom",
                    type = type,
                    user_patho = user_patho,
                    user_benign = user_benign,
                    ntree = ntree,
                    max.depth = max.depth,
                    rule.filter = 10,
                    bootstrap.rounds = lasso,
                    rule.extract.cores = cores,
                    kappa.cores = ifelse(cores < 64 & cores >= 16, 4, ifelse(cores >= 64, 10, 2)),
                    lasso.cores = ifelse(cores < 64 & cores >= 16, 5, ifelse(cores >= 64, 10, 2)))

# Save the results
saveRDS(model, paste0(scratch, "/VARPP-RuleFit_results.rds"))
varpp_report(model, paste0(scratch,"/VARPP_model"))
#=============================================================================================================================

patient_file_path  <- paste0('/mnt/', patient_vcf)
output_file_path   <- paste0('/mnt/scratch/')

#=============================================================================================================================
# Below are the fixed paths and the code

predict_data_expr <- ifelse(type %in% 'gtex', "GTEx", ifelse(type %in% 'hcl', 'HCL', 'custom'))

predict_data_expr <- ifelse(opt$type %in% 'gtex', "GTEx", ifelse(opt$type %in% 'hcl', 'HCL', 'custom'))

# Set the path for the patient file and results file
genome_file_path  <- paste0("/mnt/genome-data/", predict_data_expr) # this should be fixed either gtex or hcl (currently only gtex)
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
patient_variants <- read.vcfR(paste0('mnt/',patient_vcf))

# Bring the file into the right format (bedfile) and make sure the 'chr' is appended to the chromosome (vcf files not including this?)
as_tibble(getFIX(patient_variants)) %>%
  mutate(CHROM = paste0("chr", CHROM)) %>%
  mutate(Start = as.numeric(POS) - 1) %>%
  mutate(End = as.numeric(POS)) %>%
  select(CHROM, Start,End) -> dat

# Save the vcf so we can use bedtools
write_delim(dat, paste0(scratch,"/patient.bed"), delim="\t", col_names=F)

#================================================================================================================================
## # Run bedtools sort and merge the patient variants with the data necessary for the prediction: GTEx or HCL
system(paste0("/root/miniconda3/bin/bedtools sort -i ", scratch, "patient.bed > ", scratch, "patient_sorted.bed;
cd ", genome_file_path," ;
ls *.bed.gz | parallel '/root/miniconda3/bin/bedtools intersect -a ", scratch, "patient_sorted.bed -b {} -wb -sorted > ", scratch, "{.}_complete'"), wait=TRUE)

# Concatenate all the individual files to one result file
system(paste0("cat ", header_file, "header.txt > ", scratch,"mergelist.txt;
for i in $(seq 1 22) X Y; do cat ", scratch,"SORTED_FINAL_chr$i.bed_complete | cut -f1,2,3,8,9 --complement >>", scratch, "mergelist.txt; done
"), wait=TRUE)

# Remove the unnecessary files
#system(paste0('rm ',scratch, '*.bed_complete ', scratch, '*.bed'))

## # Read in the results file
predict <- read_delim(paste0(scratch, "mergelist.txt"), delim = "\t")
predict <- predict %>%
  mutate(CADD_raw_rankscore = CADD_2)


# Clean up all the unnecessary files and only retain the patient vcf for prediction
system(paste0("rm ", scratch, "*.bed_complete ", scratch, "patient.bed ", scratch, "patient_sorted.bed ", scratch, "mergelist.txt" ))

#================================================================================================================================
# Predict


PREDICTIONS <- predict.varppRule(predict,
                                 model,  
                                 predict="probability")

#================================================================================================================================
# Print out header of results
cat("Predictions are:\n")
head(PREDICTIONS)

# Write the results
write_csv(PREDICTIONS, paste0(scratch,"predictions.csv"))
