# R script will perform 2SMR for all exposures
# Subset exposures as necessary upstream of analysis if required

library(coloc)
library(ieugwasr)
library(genetics.binaRies)
devtools::install_github("explodecomputer/genetics.binaRies")
genetics.binaRies::get_plink_binary()
library(stringr)
library(dplyr)
library(data.table)
library(TwoSampleMR)
library(data.table)

calculate_se <- function(beta, pval) {
 
 if (length(beta) != length(pval)) {
  stop("beta and pval must have the same length")
 }
 
 z <- ifelse(
  pval < 1e-10,
  sqrt(-2 * log(pval / 2)),       
  qnorm(1 - pval / 2)            
 )
 
 z[is.infinite(z) | z == 0] <- NA
 
 se <- abs(beta / z)
 return(se)
}

# define outcome phenotype
outcome_phenotype = "" # add phenotype ID
outcome_data_dir = "/outcome/data/dir"

# define output_dir ####
output_dir = "/output/directory"

# Load and format exposure ####
exposure_filepath = "/path/to/exposure/data"
exposure_data = fread(exposure_filepath)
exposure_data_formatted = as.data.frame(matrix(nrow = nrow(exposure_data),
                                               ncol = 9))
colnames(exposure_data_formatted) = c("SNP",
                                      "beta.exposure",
                                      "se.exposure",
                                      "effect_allele.exposure",
                                      "other_allele.exposure",
                                      "eaf.exposure",
                                      "pval.exposure",
                                      "exposure",
                                      "id.exposure")

exposure_data_formatted$SNP = exposure_data$SNP
exposure_data_formatted$beta.exposure = exposure_data$beta
exposure_data_formatted$se.exposure = exposure_data$se
exposure_data_formatted$effect_allele.exposure = exposure_data$effect_allele
exposure_data_formatted$other_allele.exposure = exposure_data$other_allele
exposure_data_formatted$eaf.exposure = exposure_data$eaf
exposure_data_formatted$pval.exposure = exposure_data$pval
exposure_data_formatted$exposure = exposure_data$id.exposure
exposure_data_formatted$id.exposure = exposure_data_formatted$exposure
exposure_data$se.exposure = abs(calculate_se(exposure_data_formatted$beta.exposure, exposure_data_formatted$pval.exposure))
exposure_data_formatted$pval.exposure = as.numeric(exposure_data_formatted$pval.exposure )

# Load outcome data
outcome_data = fread("path/to/outcome/data")
outcome_data_formatted = as.data.frame(matrix(nrow = nrow(outcome_data),
                                              ncol = 9))
colnames(outcome_data_formatted) = c("SNP",
                                     "beta.outcome",
                                     "se.outcome",
                                     "effect_allele.outcome",
                                     "other_allele.outcome",
                                     "eaf.outcome",
                                     "pval.outcome",
                                     "outcome",
                                     "id.outcome")

outcome_data_formatted$SNP = outcome_data$SNP
outcome_data_formatted$beta.outcome = outcome_data$beta
outcome_data_formatted$effect_allele.outcome = outcome_data$A1
outcome_data_formatted$other_allele.outcome = outcome_data$A2
outcome_data_formatted$eaf.outcome = outcome_data$eaf
outcome_data_formatted$pval.outcome = outcome_data$pval
outcome_data_formatted$pval.outcome = as.numeric(outcome_data_formatted$pval.outcome)
outcome_data_formatted$outcome = sub(".*/", "", outcome_data$source_filepath[1])
outcome_data_formatted$id.outcome = outcome_data_formatted$outcome

outcome_data_formatted$se.outcome = calculate_se(outcome_data_formatted$beta.outcome,
                                                 outcome_data_formatted$pval.outcome)

# Define exposure QTLs ####
exposure_QTLs = unique(exposure_data_formatted$exposure)
#exposure_QTLs = exposure_data_formatted_h4snpQTLs

# Create output table ####
summary_blank_cols = c("ind", 
                       "timestamp",
                       "exposure",
                       "coloc_nsnps",
                       "mr_ivw_nsnp",
                       "mr_ivw_beta",
                       "mr_ivw_se",
                       "mr_ivw_pval",
                       "exposure_file",
                       "outcome_file")

summary_blank = as.data.frame(matrix(ncol = length(summary_blank_cols),
                                     nrow = length(exposure_QTLs)))
colnames(summary_blank) = summary_blank_cols
summary_blank$ind = 1:nrow(summary_blank)
summary_blank$exposure = exposure_QTLs
summary_blank$exposure_file[1] = sub(".*/", "", exposure_filepath)
summary_blank$outcome_file[1] = sub(".*/", "", outcome_filepath)
str(summary_blank)

output_table = summary_blank

# define output timestamp and filename
timestamp = format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
qtl_dataset = sub(".*/", "", exposure_filepath)
qtl_dataset = sub(".txt*", "", qtl_dataset)
qtl_dataset

out_dataset = sub(".*/", "", outcome_filepath)
out_dataset = sub(".txt*", "", out_dataset)
out_dataset

output_filename = paste("2SMR_",
                        qtl_dataset,
                        "_",
                        "vs",
                        "_",
                        out_dataset,
                        "_",
                        timestamp,
                        ".txt",
                        sep = "")

# Iterate over exposures ####
ind = 1
while(ind < (length(exposure_QTLs)+ 1 )){
 
 row_timestamp = format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
 
 output_table$timestamp[ind] = row_timestamp
 
 print(ind)
 
 
 
 tryCatch({
  
  ## Restrict exposure data to exposure QTL ####
  exposure_data_formatted_QTLrestr = exposure_data_formatted[exposure_data_formatted$exposure == exposure_QTLs[ind],]
  
  outcome_data_formatted_QTLrestr = outcome_data_formatted[outcome_data_formatted$SNP %in% exposure_data_formatted_QTLrestr$SNP,]
  
  mr_data = harmonise_data(exposure_data_formatted_QTLrestr,
                           outcome_data_formatted_QTLrestr,
                           action = 1)
  str(mr_data)
  
  mr_res = mr(mr_data)
  
  if ("Inverse variance weighted" %in% mr_res$method) {
   mr_subset = subset(mr_res, mr_res$method == "Inverse variance weighted")
   output_table$mr_ivw_nsnp[ind] = mr_subset$nsnp
   output_table$mr_ivw_beta[ind] = mr_subset$b
   output_table$mr_ivw_se[ind] = mr_subset$se
   output_table$mr_ivw_pval[ind] = mr_subset$pval
  }
  
  write.table(output_table,
              file = paste(output_dir, output_filename, sep = "/"),
              quote = FALSE,
              col.names = TRUE,
              row.names = FALSE)
  
  
 }, error = function(e) {
  message("An error occurred: ", e$message)
 })
 
 ind = ind+1
 
}