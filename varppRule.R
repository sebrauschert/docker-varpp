#!/usr/bin Rscript

# Script to run VARPP-RuleFit 
# 

# Load the necessary libraries
library(optparse)
#=============================================================================================================================

# OPTPARSE options
option_list = list(
  make_option(c("-g", "--gene_list"),
              type="character",
              default=NULL,
              help="A list of potentially pathogenic genes associated with the patient.",
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
              help="Number of LASSO bootstrap rounds.",
              metavar="character"),
  make_option(c("-m", "--max_depth"),
              type="character",
              default=3,
              help="Number of LASSO bootstrap rounds.",
              metavar="character")
);

opt_parser <- OptionParser(option_list = option_list,
                           description = "\n VARPP-RuleFit.",
                           epilogue = "Example:\n\n  ./varppRule.R -g <gene list> -t <type: GTEx, HCL, custom> -p <User provided pathogenic data> -b <User provided benign data> -o <output/path> -c <Maximum number of cores available> -n <Number of trees, defaults to 2000> -l <Number of LASSO bootstrap rounds> -m <Tree depth>\n\n");

opt <- parse_args(opt_parser)

library(varppRule)
library(tidyverse)

genes       <- str_trim(unlist(str_split(opt$gene_list, ",")))
type        <- opt$type
user_patho  <- opt$user_pathogenic
user_benign <- opt$user_benign
output_path <- opt$output_path
cores       <- as.numeric(opt$cores)
ntree       <- as.numeric(opt$ntree)
max.depth   <- as.numeric(opt$max_depth)
lasso       <- as.numeric(opt$lasso)

#genes <-  "ATL3,COL3A1,CTSC,CTSK,DDR2"
# genes <- c('ATL3', 'COL3A1', 'CTSC', 'CTSK', 'DDR2')
# type  <- 'gtex'
# ntree <- 200
# cores <- 64



# Run the model
results <- rule_fit(HPO_genes = genes,
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
saveRDS(results, paste0(output_path, "/VARPP-RuleFit_results.rds"))