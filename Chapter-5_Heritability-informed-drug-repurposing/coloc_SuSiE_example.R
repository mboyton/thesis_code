library(LAVA)
library(coloc)
library(ieugwasr)
#devtools::install_github("explodecomputer/genetics.binaRies")
library(genetics.binaRies)
genetics.binaRies::get_plink_binary()
library(stringr)
library(dplyr)
library(data.table)
#library(TwoSampleMR)

output_directory = "output/directory"

# Define univ p_thresh ####
bivar_p_thresh = 0.05
bivariate_id = "id_for_output" # change this as appropriate

# Bivariate LAVA results directory
bivar_results_dir = "/path/to/bivariate/results"
bivar_results_dir_contents = list.files(bivar_results_dir)

# Define and load phenotypes

phenotype_1 = "" # amend as appropriate
phenotype_2 = "" # amend as appropriate

phenotype_1_raw_data = fread("path/to/phenotype_1/gwas/data", fill = TRUE)

phenotype_2_raw_data = fread("path/to/phenotype_2/gwas/data", fill = TRUE)

# load reference SNP dataset for coloc
ref_snps_all = read.table("path/to/EUR.bim",
                          header=FALSE) # 1000G TwoSampleMR reference
colnames(ref_snps_all) = c("chr", "rsid", ".", "..", "A1", "A2")

# define LAVA dependencies
input_info_file = "path/to/input/info"
sample_overlap_file = "path/to/sample/overlap/file"
ref_prefix = "path/to/g1000_eur"
loci = "path/to/loci_2495.txt"

input = process.input(input.info.file = input_info_file,           # input info file
                      sample.overlap.file = sample_overlap_file,   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix = ref_prefix,                    # reference genotype data prefix
                      phenos = c(phenotype_1, phenotype_2))

loci = read.loci(loci)

# load bivariate results
bivar_results = fread("path/to/bivariate/results")

# subset bivar_results according to p_thresh
bivar_results = subset(bivar_results,
                       bivar_results$p < bivar_p_thresh)

n_bivar_results = nrow(bivar_results)
print(paste(n_bivar_results, "significant LAVA loci between phenotype pair"))

# loop through significant loci
loci_for_susie = bivar_results$LOC

# generate output log file
output_log_colnames = c("locus_i",
                        "timestamp",
                        "locus",
                        "chr",
                        "start",
                        "stop",
                        "locus_snps",
                        "p1_snps",
                        "p2_snps",
                        "p1p2_snps",
                        "ref_snps",
                        "p1_ref_mismatch",
                        "p1_ref_req_align",
                        "p2_ref_mismatch",
                        "p2_ref_req_align",
                        "p1_beta_zero",
                        "p2_beta_zero",
                        "p1_betase_zero",
                        "p2_betase_zero",
                        "p1_snps_not_in_ldmatrix",
                        "p2_snps_not_in_ldmatrix",
                        "ldmatrix_nsnps",
                        "phenotype_1",
                        "phenotype_2"
)

output_log = as.data.frame(matrix(ncol = length(output_log_colnames),
                                  nrow = length(loci_for_susie)))
colnames(output_log) = output_log_colnames
output_log$locus_i = 1:length(loci_for_susie)

# define timestamp for coloc SuSiE output
current_time = Sys.time()
timestamp = format(current_time, "%Y_%m_%d_%H_%M_%S")

susie.res.list = list()

# iterate across loci
locus_i = 1

while(locus_i < length(loci_for_susie)+1){
 
 
 locus_index = loci_for_susie[locus_i]
 locus = process.locus(loci[locus_index ,], input)
 
 # update the log
 output_log[locus_i,]$locus_i = locus_i
 output_log[locus_i,]$locus = locus_index
 output_log[locus_i,]$chr = loci[locus_index ,]$CHR
 output_log[locus_i,]$start = loci[locus_index ,]$START
 output_log[locus_i,]$stop = loci[locus_index ,]$STOP
 output_log[locus_i,]$locus_snps = length(locus$snps)
 
 # subset raw data files for SNPs in locus
 phenotype_1_raw_data_snps = phenotype_1_raw_data$SNP %in% locus$snps
 phenotype_1_raw_data_subset = phenotype_1_raw_data[phenotype_1_raw_data_snps,]
 
 print(paste("Phenotype 1:", nrow(phenotype_1_raw_data_subset), "SNPs available in summary data"))
 
 phenotype_2_raw_data_snps = phenotype_2_raw_data$SNP %in% locus$snps
 phenotype_2_raw_data_subset = phenotype_2_raw_data[phenotype_2_raw_data_snps,]
 
 print(paste("Phenotype 2:", nrow(phenotype_2_raw_data_subset), "SNPs available in summary data"))
 
 # update the log
 output_log[locus_i,]$p1_snps = nrow(phenotype_1_raw_data_subset)
 output_log[locus_i,]$p2_snps = nrow(phenotype_2_raw_data_subset)
 
 # compile vector of SNPs required
 pheno_snps_loci_all = unique(phenotype_1_raw_data_subset$SNP, phenotype_2_raw_data_subset$SNP)
 
 print(paste(length(pheno_snps_loci_all), " SNPs shared between phenotypes at locus"))
 
 # update the log
 output_log[locus_i,]$p1p2_snps = length(pheno_snps_loci_all)
 
 # subset ref_snps to include only SNPs required
 snps_in_ref = (ref_snps_all$rsid %in% pheno_snps_loci_all)
 ref_snps = ref_snps_all[(snps_in_ref),]
 print(paste(nrow(ref_snps), " SNPs available in reference dataset"))
 
 # update the log
 output_log[locus_i,]$ref_snps = nrow(ref_snps)
 
 if(nrow(ref_snps) == 0){
  locus_i = locus_i + 1
 }
 
 if(nrow(ref_snps) > 0){
  # subset pheno datasets to include only SNPs which are present in ref_snps 
  phenotype_1_raw_data_subset = phenotype_1_raw_data_subset[phenotype_1_raw_data_subset$SNP %in% ref_snps$rsid, ]
  phenotype_2_raw_data_subset = phenotype_2_raw_data_subset[phenotype_2_raw_data_subset$SNP %in% ref_snps$rsid, ]
  
  
  # order the phenotype_raw_data_subsets according to ref_snps$rsid
  phenotype_1_raw_data_subset = phenotype_1_raw_data_subset[match(ref_snps$rsid, phenotype_1_raw_data_subset$SNP),]
  phenotype_2_raw_data_subset = phenotype_2_raw_data_subset[match(ref_snps$rsid, phenotype_2_raw_data_subset$SNP),]
  
  head(phenotype_1_raw_data_subset)
  head(phenotype_2_raw_data_subset)
  
  
  summary(phenotype_1_raw_data_subset$SNP == ref_snps$rsid) # qc
  summary(phenotype_2_raw_data_subset$SNP == ref_snps$rsid) # qc
  
  # --- Harmonise PHENOTYPE_1 to the reference (vectorised) ---
  # Assumes: phenotype_1_raw_data_subset has columns SNP, A1, A2, beta
  #          ref_snps has columns rsid (SNP id), A1, A2
  
  # 1) Align by SNP to the reference
  m <- match(phenotype_1_raw_data_subset$SNP, ref_snps$rsid)
  keep_in_ref <- !is.na(m)
  
  p1  <- phenotype_1_raw_data_subset[keep_in_ref, , drop = FALSE]
  ref <- ref_snps[m[keep_in_ref], , drop = FALSE]
  
  # 2) Classify relationships vs reference
  is_match <- (p1$A1 == ref$A1 & p1$A2 == ref$A2)
  is_flip  <- (p1$A1 == ref$A2 & p1$A2 == ref$A1)
  is_err   <- !(is_match | is_flip)  # present in ref but not resolvable
  
  # 3) Flip alleles / beta where needed
  if (any(is_flip)) {
   tmpA1        <- p1$A1[is_flip]
   p1$A1[is_flip]   <- p1$A2[is_flip]
   p1$A2[is_flip]   <- tmpA1
   p1$beta[is_flip] <- -p1$beta[is_flip]
  }
  
  # 4) Drop irreconcilable rows
  removed <- sum(is_err)
  p1_ok   <- p1[!is_err, , drop = FALSE]
  
  # 5) Write back: replace kept SNPs in the original table and drop the rest
  # (overwrites only harmonised rows; removes SNPs not kept)
  phenotype_1_raw_data_subset <- p1_ok
  
  # 6) Log + message
  output_log[locus_i, ]$p1_ref_mismatch  <- removed
  output_log[locus_i, ]$p1_ref_req_align <- nrow(p1_ok)
  print(paste(removed, "SNPs removed that do not match the reference dataset"))
  
  # --- Harmonise PHENOTYPE_2 to the reference (vectorised) ---
  # Assumes: phenotype_2_raw_data_subset has columns SNP, A1, A2, beta
  #          ref_snps has columns rsid, A1, A2
  
  # 1) Align by SNP to the reference
  m2 <- match(phenotype_2_raw_data_subset$SNP, ref_snps$rsid)
  keep_in_ref2 <- !is.na(m2)
  
  p2  <- phenotype_2_raw_data_subset[keep_in_ref2, , drop = FALSE]
  ref <- ref_snps[m2[keep_in_ref2], , drop = FALSE]
  
  # 2) Classify vs reference
  is_match2 <- (p2$A1 == ref$A1 & p2$A2 == ref$A2)
  is_flip2  <- (p2$A1 == ref$A2 & p2$A2 == ref$A1)
  is_err2   <- !(is_match2 | is_flip2)
  
  # 3) Flip alleles / beta where needed
  if (any(is_flip2)) {
   tmpA1           <- p2$A1[is_flip2]
   p2$A1[is_flip2] <- p2$A2[is_flip2]
   p2$A2[is_flip2] <- tmpA1
   p2$beta[is_flip2] <- -p2$beta[is_flip2]
  }
  
  # 4) Drop irreconcilable rows
  removed2 <- sum(is_err2)
  p2_ok    <- p2[!is_err2, , drop = FALSE]
  
  # 5) Write back: keep only harmonised rows
  phenotype_2_raw_data_subset <- p2_ok
  
  # 6) Log + message
  output_log[locus_i, ]$p2_ref_mismatch  <- removed2
  output_log[locus_i, ]$p2_ref_req_align <- nrow(p2_ok)
  print(paste(removed2, "SNPs removed that do not match the reference dataset"))
  
  # ---- remove SNPs with beta == 0 in either phenotype, then re-sync ----
  # find zero-beta SNPs
  pheno_1_snps_zero <- phenotype_1_raw_data_subset$SNP[phenotype_1_raw_data_subset$beta == 0]
  pheno_2_snps_zero <- phenotype_2_raw_data_subset$SNP[phenotype_2_raw_data_subset$beta == 0]
  
  # log counts
  output_log[locus_i, ]$p1_beta_zero <- length(pheno_1_snps_zero)
  output_log[locus_i, ]$p2_beta_zero <- length(pheno_2_snps_zero)
  
  # union of zero-beta SNPs across both traits
  snps_to_remove_before_coloc <- unique(c(pheno_1_snps_zero, pheno_2_snps_zero))
  print(paste(length(snps_to_remove_before_coloc), " SNPs removed with beta = 0"))
  
  # drop from both phenotypes (use SNP, not snp)
  if (length(snps_to_remove_before_coloc) > 0) {
   keep1 <- !(phenotype_1_raw_data_subset$SNP %in% snps_to_remove_before_coloc)
   keep2 <- !(phenotype_2_raw_data_subset$SNP %in% snps_to_remove_before_coloc)
   phenotype_1_raw_data_subset <- phenotype_1_raw_data_subset[keep1, , drop = FALSE]
   phenotype_2_raw_data_subset <- phenotype_2_raw_data_subset[keep2, , drop = FALSE]
  }
  
  # re-sync by SNP (ensure identical set and order in both tables)
  snps_common <- intersect(phenotype_1_raw_data_subset$SNP, phenotype_2_raw_data_subset$SNP)
  phenotype_1_raw_data_subset <- phenotype_1_raw_data_subset[phenotype_1_raw_data_subset$SNP %in% snps_common, , drop = FALSE]
  phenotype_2_raw_data_subset <- phenotype_2_raw_data_subset[phenotype_2_raw_data_subset$SNP %in% snps_common, , drop = FALSE]
  
  # order phenotype_2 to match phenotype_1 SNP order
  o2 <- match(phenotype_1_raw_data_subset$SNP, phenotype_2_raw_data_subset$SNP)
  phenotype_2_raw_data_subset <- phenotype_2_raw_data_subset[o2, , drop = FALSE]
  
  # QC (SNPs and alleles should match 1:1)
  summary(phenotype_1_raw_data_subset$SNP == phenotype_2_raw_data_subset$SNP)
  summary(phenotype_1_raw_data_subset$A1  == phenotype_2_raw_data_subset$A1)
  summary(phenotype_1_raw_data_subset$A2  == phenotype_2_raw_data_subset$A2)
  
  # calculate beta_se for both phenotypes
  calculate_se <- function(beta, pval) {
   # Calculate Z-score from p-value
   Z <- qnorm(1 - pval / 2)
   
   # Calculate standard error
   SE <- beta / Z
   
   return(SE)
  }
  
  # ---- calculate SE and drop zero/Inf ----
  
  # 1) calculate beta_se
  phenotype_1_raw_data_subset$beta_se <- calculate_se(
   beta = phenotype_1_raw_data_subset$beta,
   pval = as.numeric(phenotype_1_raw_data_subset$pval)
  )
  
  phenotype_2_raw_data_subset$beta_se <- calculate_se(
   beta = phenotype_2_raw_data_subset$beta,
   pval = as.numeric(phenotype_2_raw_data_subset$pval)
  )
  
  # 2) identify bad SNPs
  p1_bad <- phenotype_1_raw_data_subset$SNP[
   phenotype_1_raw_data_subset$beta_se == 0 |
    is.infinite(phenotype_1_raw_data_subset$beta_se)
  ]
  
  p2_bad <- phenotype_2_raw_data_subset$SNP[
   phenotype_2_raw_data_subset$beta_se == 0 |
    is.infinite(phenotype_2_raw_data_subset$beta_se)
  ]
  
  # 3) update the log
  output_log[locus_i, ]$p1_betase_zero <- length(p1_bad)
  output_log[locus_i, ]$p2_betase_zero <- length(p2_bad)
  
  # 4) remove bad SNPs from both phenotypes
  bad_all <- unique(c(p1_bad, p2_bad))
  phenotype_1_raw_data_subset <- phenotype_1_raw_data_subset[!(phenotype_1_raw_data_subset$SNP %in% bad_all), ]
  phenotype_2_raw_data_subset <- phenotype_2_raw_data_subset[!(phenotype_2_raw_data_subset$SNP %in% bad_all), ]
  
  # 5) re-sync SNP sets between phenotypes
  snps_common <- intersect(phenotype_1_raw_data_subset$SNP, phenotype_2_raw_data_subset$SNP)
  phenotype_1_raw_data_subset <- phenotype_1_raw_data_subset[phenotype_1_raw_data_subset$SNP %in% snps_common, ]
  phenotype_2_raw_data_subset <- phenotype_2_raw_data_subset[phenotype_2_raw_data_subset$SNP %in% snps_common, ]
  
  # 6) QC alignment
  summary(phenotype_1_raw_data_subset$SNP == phenotype_2_raw_data_subset$SNP)
  summary(phenotype_1_raw_data_subset$A1  == phenotype_2_raw_data_subset$A1)
  summary(phenotype_1_raw_data_subset$A2  == phenotype_2_raw_data_subset$A2)
  
  # add LD information for coloc susie ####
  snps_for_ld_matrix = phenotype_1_raw_data_subset$SNP
  snps_ld_matrix = ld_matrix(
   snps_for_ld_matrix,
   plink_bin = genetics.binaRies::get_plink_binary(),
   bfile = "path/to/bfile/EUR"
  )
  
  
  snps_in_ld_matrix = unique(c(colnames(snps_ld_matrix),
                               rownames(snps_ld_matrix)))
  
  snps_in_ld_matrix.rm_alleles = str_extract(snps_in_ld_matrix, "[^_]+")
  
  # subset phenotype datasets to include data with LD reference only
  phenotype_1_raw_data_subset_forSUSIE = phenotype_1_raw_data_subset[(phenotype_1_raw_data_subset$SNP %in% snps_in_ld_matrix.rm_alleles),]
  phenotype_2_raw_data_subset_forSUSIE = phenotype_2_raw_data_subset[(phenotype_2_raw_data_subset$SNP %in% snps_in_ld_matrix.rm_alleles),]
  
  ld_ref_qc = as.data.frame(matrix(nrow = length(snps_in_ld_matrix)))
  
  ld_ref_qc$snp = snps_in_ld_matrix.rm_alleles
  ld_ref_qc$ld_matrix_A1 <- sub("^.*_(.)_(.)$", "\\1", snps_in_ld_matrix)  # Extract the first allele
  ld_ref_qc$ld_matrix_A2 <- sub("^.*_(.)_(.)$", "\\2", snps_in_ld_matrix)  # Extract the second allele
  
  ld_ref_qc$phenos_A1 = phenotype_1_raw_data_subset_forSUSIE$A1
  ld_ref_qc$phenos_A2 = phenotype_1_raw_data_subset_forSUSIE$A2
  
  ld_ref_qc$A1_check = as.numeric(ld_ref_qc$ld_matrix_A1 == ld_ref_qc$phenos_A1)
  ld_ref_qc$A2_check = as.numeric(ld_ref_qc$ld_matrix_A2 == ld_ref_qc$phenos_A2)
  
  ld_ref_qc <- ld_ref_qc %>%
   mutate(
    alignment_check = case_when(
     ld_matrix_A1 == phenos_A1 & ld_matrix_A2 == phenos_A2 ~ "match",
     ld_matrix_A1 == phenos_A2 & ld_matrix_A2 == phenos_A1 ~ "flip",
     TRUE ~ "error"
    )
   )
  
  # remove error SNPs from ld reference panel and phenotype datasets
  snps_failed_qc_ref = ld_ref_qc$alignment_check == "error"
  phenotype_1_raw_data_subset_forSUSIE = phenotype_1_raw_data_subset_forSUSIE[!snps_failed_qc_ref,]
  phenotype_2_raw_data_subset_forSUSIE = phenotype_2_raw_data_subset_forSUSIE[!snps_failed_qc_ref,]
  
  snps_for_coloc_ld_matrix = ld_ref_qc[!snps_failed_qc_ref,]
  snps_for_coloc_ld_matrix$multiplier = snps_for_coloc_ld_matrix$alignment_check
  snps_for_coloc_ld_matrix$multiplier[snps_for_coloc_ld_matrix$multiplier == "match"] = "1"
  snps_for_coloc_ld_matrix$multiplier[snps_for_coloc_ld_matrix$multiplier == "flip"] = "-1"
  snps_for_coloc_ld_matrix$multiplier = as.numeric(snps_for_coloc_ld_matrix$multiplier)
  
  snp_multiplier = snps_for_coloc_ld_matrix[,c("snp", "multiplier")]
  
  snps_ld_matrix = snps_ld_matrix[!snps_failed_qc_ref,!snps_failed_qc_ref]
  
  # update the log
  output_log[locus_i,]$p1_snps_not_in_ldmatrix = nrow(phenotype_1_raw_data_subset_forSUSIE[snps_failed_qc_ref,])
  output_log[locus_i,]$p2_snps_not_in_ldmatrix = nrow(phenotype_1_raw_data_subset_forSUSIE[snps_failed_qc_ref,])
  
  
  # Convert dataframe to matrix
  ld_multiplier <- matrix(0, nrow = nrow(snp_multiplier), ncol = nrow(snp_multiplier))
  rownames(ld_multiplier) <- snp_multiplier$snp
  colnames(ld_multiplier) <- snp_multiplier$snp
  
  # Fill in the matrix with products
  for (i in 1:nrow(ld_multiplier)) {
   for (j in 1:ncol(ld_multiplier)) {
    ld_multiplier[i, j] <- snp_multiplier$multiplier[i] * snp_multiplier$multiplier[j]
   }
  }
  
  # convert the snps_ld_matrix
  snps_ld_matrix = snps_ld_matrix * ld_multiplier
  colnames(snps_ld_matrix) = str_extract(colnames(snps_ld_matrix), "[^_]+")
  rownames(snps_ld_matrix) = str_extract(rownames(snps_ld_matrix), "[^_]+")
  
  # update the log
  output_log[locus_i,]$ldmatrix_nsnps = nrow(snps_ld_matrix)
  
  # construct SuSiE data objects
  
  ## retrieve the study samples
  input_info_file_df = read.table(input_info_file,
                                  header=TRUE)
  
  phenotype_1_input_info = subset(input_info_file_df,
                                  input_info_file_df$phenotype == phenotype_1)
  phenotype_2_input_info = subset(input_info_file_df,
                                  input_info_file_df$phenotype == phenotype_2)
  
  # create SuSiE objects
  data_phenotype_1_SUSIE = list()
  data_phenotype_1_SUSIE$beta = phenotype_1_raw_data_subset_forSUSIE$beta
  data_phenotype_1_SUSIE$varbeta = phenotype_1_raw_data_subset_forSUSIE$beta_se^2
  data_phenotype_1_SUSIE$snp = phenotype_1_raw_data_subset_forSUSIE$SNP
  data_phenotype_1_SUSIE$position = c(1:nrow(phenotype_1_raw_data_subset_forSUSIE))
  data_phenotype_1_SUSIE$type = "cc"
  data_phenotype_1_SUSIE$LD = snps_ld_matrix
  data_phenotype_1_SUSIE$N = phenotype_1_input_info$cases + phenotype_1_input_info$controls
  
  data_phenotype_2_SUSIE = list()
  data_phenotype_2_SUSIE$beta = phenotype_2_raw_data_subset_forSUSIE$beta
  data_phenotype_2_SUSIE$varbeta = phenotype_2_raw_data_subset_forSUSIE$beta_se^2
  data_phenotype_2_SUSIE$snp = phenotype_2_raw_data_subset_forSUSIE$SNP
  data_phenotype_2_SUSIE$position = c(1:nrow(phenotype_2_raw_data_subset_forSUSIE))
  data_phenotype_2_SUSIE$type = "cc"
  data_phenotype_2_SUSIE$LD = snps_ld_matrix
  data_phenotype_2_SUSIE$N = phenotype_2_input_info$cases + phenotype_2_input_info$controls
  
  t <- runsusie(data_phenotype_1_SUSIE)
  tt <- runsusie(data_phenotype_2_SUSIE)
  
  
  susie_output_filename = paste(bivariate_id,
                                "_susie_res_",
                                timestamp,
                                ".rda", sep = "")
  
  susie.res <- coloc.susie(t, tt)
  list_entry_name = as.character(locus_i)
  susie.res.list[[list_entry_name]] = susie.res
  saveRDS(susie.res.list, file = paste(output_directory, "/", susie_output_filename, sep = ""))
  
  # update the log
  output_log[locus_i,]$timestamp = format(current_time, "%Y_%m_%d_%H_%M_%S")
  output_log[locus_i,]$phenotype_1 = phenotype_1
  output_log[locus_i,]$phenotype_2 = phenotype_2
  
  write.table(output_log,
              file = paste(output_directory,
                           "/",
                           bivariate_id,
                           "_output_log_",
                           timestamp,
                           ".txt", sep = ""),
              quote = FALSE,
              row.names = FALSE)
  
  locus_i <- locus_i + 1
 }
}
