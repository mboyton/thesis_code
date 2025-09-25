# Load Libraries
library("TwoSampleMR")
library("readxl")
library("SciViews")
library(ggplot2)
library(dplyr)
library(patchwork)
library(stringr)

# 1. Load and format exposure data ####
exp_data_filepath <- "/path/to/exposure/data.txt"
exposure_data <- read.delim(exp_data_filepath)

## define colnames
colnames(exposure_data) <- c("phenotype",
                             "SNP",
                             "chr_name",
                             "chrom_start",
                             "effect_allele.exposure",
                             "other_allele.exposure",
                             "beta.exposure",
                             "se.exposure",
                             "pval")
exposure_data$id.exposure <- exposure_data$phenotype
exposure_data$exposure <- exposure_data$phenotype
exposure_data$eaf.exposure <- exposure_data$eaf
exposure_data$SNP_id <- str_c(exposure_data$chr_name, ":",
                              exposure_data$chrom_start, ":",
                              exposure_data$effect_allele.exposure, ":",
                              exposure_data$other_allele.exposure)

# 2. Load and format outcome data ####
load("/path/to/outcome/data.Robj", verbose=TRUE)


# 3. Perform TwoSampleMR ####
## Create output dataframe mr_output_df ####
mr_output_df <- as.data.frame(matrix(ncol=9,nrow=0))
colnames(mr_output_df) <- c("id.exposure",
                            "id.outcome",
                            "outcome",
                            "exposure",
                            "method",
                            "nsnp",
                            "b",
                            "se",
                            "pval")

exposures <- unique(exposure_data$phenotype)

## Iterate over exposures ####
i <- 1
while(i<length(exposures)+1){
 
 print(paste("Analysing exposure", i, "of", length(exposures), sep = " "))
 exp <- subset(exposure_data, exposure_data$phenotype == exposures[i])
 snps_include <- out_dat$SNP %in% exp$SNP
 out <- out_dat[snps_include,]
 
 if(nrow(exp) == 0){
  i <- i+1
 }else if(nrow(exp) > 1){
  
  dat <- harmonise_data(exp, out, action=2)
  output <- mr(dat)
  mr_output_df <- rbind(mr_output_df, output)
  i <- i+1
 }
}

# 4. Filter for IVW results and apply FDR adjustment ####
mr_output_df$gene <- gsub(".*\\_", "", mr_output_df$id.exposure)

mr_output_IVW_df <- subset(mr_output_df,
                           mr_output_df$method == "Inverse variance weighted")
mr_output_IVW_df <- mr_output_IVW_df[order(mr_output_IVW_df$pval),]
mr_output_IVW_df$pval.fdr <- p.adjust(mr_output_IVW_df$pval, "fdr")
mr_output_IVW_df_sig <- subset(mr_output_IVW_df, pval.fdr < 0.05)

# 5. Define comorbid outcomes ####
# comorbidities opengwas IDs:
depression = "ieu-b-102"
anxiety = "finn-b-F5_GAD"
hypertension = "ukb-b-12493"
hypercholesterolemia = "finn-b-E4_HYPERCHOL"
## chronic lung diseases
asthma = "ebi-a-GCST90014325"
COPD = "finn-b-J10_COPD"

comorb_outcomes <- c(depression, anxiety, hypertension, hypercholesterolemia, asthma, COPD)

outcomes <- comorb_outcomes

# retrieve data from OpenGWAS
out_data_all <- extract_outcome_data(exposure_data$SNP, outcomes = outcomes)

# 6. Perform comorbidity-focused MR ####
# add multiphenome_output_df
exposures <- mr_output_IVW_df_sig$id.exposure
multiphenome_output_df <- as.data.frame(matrix(ncol=9,nrow=0))
colnames(multiphenome_output_df) <- c("id.exposure",
                                      "id.outcome",
                                      "outcome",
                                      "exposure",
                                      "method",
                                      "nsnp",
                                      "b",
                                      "se",
                                      "pval")


phewas_generate <- function(exposure_name,
                            output_filepath){
 
 ## re-init inside the function so results don't leak across exposures
 multiphenome_output_df <- as.data.frame(matrix(ncol=9,nrow=0))
 colnames(multiphenome_output_df) <- c("id.exposure",
                                       "id.outcome",
                                       "outcome",
                                       "exposure",
                                       "method",
                                       "nsnp",
                                       "b",
                                       "se",
                                       "pval")
 
 ii <- 1
 while(ii<length(outcomes)+1){
  print(paste("Analysing against", outcomes[ii], "outcome", ii, "of", length(outcomes), sep = " "))
  
  exp_dat <- subset(exposure_data, exposure_data$phenotype == exposure_name)
  out_dat <- subset(out_data_all, out_data_all$id.outcome == outcomes[ii])
  
  snps_keep <- out_dat$SNP %in% exp_dat$SNP
  out_dat <- out_dat[snps_keep,]
  
  dat <- harmonise_data(exp_dat, out_dat, action=2)
  
  ## only run MR if at least one SNP remains after harmonisation
  if(nrow(dat) > 0){
   results <- mr(dat)
   multiphenome_output_df <- rbind(multiphenome_output_df, results)
  }
  
  ii <- ii + 1
 }
 
 output_filename <- paste(exposure_name, "PhewasOutput", sep = "_")
 output_filename <- paste(output_filename, "Robj", sep = ".")
 output_filename <- paste(output_filepath, output_filename, sep = "/")
 save(multiphenome_output_df, file = output_filename)
}

i <- 1
while(i < length(exposures)+1){
 
 
 phewas_generate(exposure_name <- exposures[i],
                 output_filepath = "/path/to/output/directory")
 i <- i+1
}
