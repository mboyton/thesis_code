library(LAVA)
library(stringr)
library(data.table)

setwd("/path/to/working/directory")

univariate_data_directory = "path/to/univariate/data"

# Define univ p_thresh
univ_p_thresh = 0.05

# Define phenotype pair ####

phenotype_1 = "" # define as appropriate
phenotype_2 = "" # define as appropriate

# create temporary sample_overlap file
sample_overlap_temp_dir = getwd()
sample_overlap_temp_id = bivariate_id
sample_overlap_temp_id = paste("analysis", sample_overlap_temp_id, sep = "_")
sample_overlap_temp_id = paste(sample_overlap_temp_id, "sampleoverlap.txt", sep = "_")

gcov_int_matrix = read.table("/path/to/matrix/of/gcov_int/values.txt")

gcov_int_matrix = as.matrix(gcov_int_matrix)

gcov_int_matrix_subset = gcov_int_matrix[rownames(gcov_int_matrix) %in% c(phenotype_1,
                                                                          phenotype_2),
                                         colnames(gcov_int_matrix) %in% c(phenotype_1,
                                                                          phenotype_2)]
# Set diagonal values to 1
diag(gcov_int_matrix_subset) = 1
gcov_int_matrix_subset
write.table(gcov_int_matrix_subset, file = sample_overlap_temp_id, row.names=TRUE,sep="\t", quote = FALSE)

# load univariate data for both phenotypes

pheno_1_univ_res_dir = paste(univariate_data_directory,
                             phenotype_1,
                             sep = "/")
pheno_1_univ_res_dir_files = list.files(pheno_1_univ_res_dir)
pheno_1_univ_res = pheno_1_univ_res_dir_files[grep("LAVA_univ", pheno_1_univ_res_dir_files)]
pheno_1_univ_res_filepath = paste(pheno_1_univ_res_dir,
                                  "/",
                                  pheno_1_univ_res,
                                  sep = "")

pheno_2_univ_res_dir = paste(univariate_data_directory,
                             phenotype_2,
                             sep = "/")
pheno_2_univ_res_dir_files = list.files(pheno_2_univ_res_dir)
pheno_2_univ_res = pheno_2_univ_res_dir_files[grep("LAVA_univ", pheno_2_univ_res_dir_files)]
pheno_2_univ_res_filepath = paste(pheno_2_univ_res_dir,
                                  pheno_2_univ_res,
                                  sep = "/")


pheno_1_univ = fread(pheno_1_univ_res_filepath)
pheno_2_univ = fread(pheno_2_univ_res_filepath)

# Load LAVA dependencies ####
# define LAVA dependencies
input_info_file = "/path/to/input/info/file.txt"
sample_overlap_file = sample_overlap_temp_id
ref_prefix = "/path/to/g1000_eur"
loci = "/path/to/loci_2495.txt"


input = process.input(input.info.file = input_info_file,           # input info file
                      sample.overlap.file = sample_overlap_file,   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix = ref_prefix,                    # reference genotype data prefix
                      phenos = c(phenotype_1, phenotype_2))

loci = read.loci(loci)

# Determine intersect of significant univariate signals ####
pheno_1_sigloci = pheno_1_univ[(pheno_1_univ$p < univ_p_thresh),]$LOC

pheno_2_sigloci = pheno_2_univ[(pheno_2_univ$p < univ_p_thresh),]$LOC

common_loci = intersect(pheno_1_sigloci, pheno_2_sigloci)

print(paste("Total common loci = ", length(common_loci), " of 2495", sep = ""))

# Create bivariate output table ####
bivar.table <- as.data.frame(matrix(ncol=13, nrow = nrow(loci)))
colnames(bivar.table) <- c("LOC","CHR","START","STOP","phen1","phen2", "rho", "rho.lower", "rho.upper", "r2", "r2.lower", "r2.upper", "p")
bivar.table$LOC = loci$LOC
bivar.table$CHR = loci$CHR
bivar.table$START = loci$START
bivar.table$STOP = loci$STOP

# Iterate LAVA through common_loci ####

i = 1
while(i<length(common_loci)+1){
 
 print(paste("Processing locus ", i, " of ", length(common_loci), sep = ""))
 
 locus = process.locus(loci[common_loci[i],], input)
 bivar.data = run.bivar(locus, param.lim = 1000)
 
 bivar.table[bivar.table$LOC == common_loci[i],]$phen1 = bivar.data$phen1
 bivar.table[bivar.table$LOC == common_loci[i],]$phen2 = bivar.data$phen2
 bivar.table[bivar.table$LOC == common_loci[i],]$rho = bivar.data$rho
 bivar.table[bivar.table$LOC == common_loci[i],]$rho.lower = bivar.data$rho.lower
 bivar.table[bivar.table$LOC == common_loci[i],]$rho.upper = bivar.data$rho.upper
 bivar.table[bivar.table$LOC == common_loci[i],]$r2 = bivar.data$r2
 bivar.table[bivar.table$LOC == common_loci[i],]$r2.lower = bivar.data$r2.lower
 bivar.table[bivar.table$LOC == common_loci[i],]$r2.upper = bivar.data$r2.upper
 bivar.table[bivar.table$LOC == common_loci[i],]$p = bivar.data$p
 
 i = i+1
 
}

timestamp = format(Sys.time(), "%Y_%m_%d_%H_%M_%S")

filename = paste(bivariate_id,
                 "_",
                 phenotype_1,
                 "_",
                 phenotype_2,
                 "_",
                 "LAVA_bivar",
                 "_",
                 timestamp,
                 ".txt",
                 sep = "")


write.table(bivar.table,
            file = filename,
            row.names=FALSE,
            sep="\t",
            quote = FALSE)

siglociinfo = bivar.table[!is.na(bivar.table$p),]
siglociinfo = nrow(siglociinfo[(siglociinfo$p < 0.05), ])
siglociinfo

print(paste("Bivariate LAVA complete"))
print(paste(siglociinfo, " loci significant at p < 0.05", sep = ""))