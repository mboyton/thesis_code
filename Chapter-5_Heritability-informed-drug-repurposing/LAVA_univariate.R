library(LAVA)
library(data.table)

# define phenotype
phenotype = "phenotype_name" # adjust as appropriate

# define LAVA dependencies
input_info_file = "/path/to/input_info_file.txt"
sample_overlap_file = NULL
ref_prefix = "/path/to/g1000_eur"
loci = "/path/to/loci_2495.txt"

# generate univariate LAVA results

input = process.input(input.info.file=input_info_file,           # input info file
                      sample.overlap.file=sample_overlap_file,   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix=ref_prefix,                    # reference genotype data prefix
                      phenos=phenotype)

loci = read.loci(loci)

chroms = unique(loci$CHR)

chr = 1 # initate loop ####
while(chr<length(chroms)+1){ # length(chroms)+1
 # while(chr<length(chroms)+1){ # length(chroms)+1
 
 print(paste("###### ANALYSING CHROMOSOME", chr, "OF", length(chroms), "#####", sep = " "))
 loci_analyse = subset(loci, loci$CHR == chr) # change this in loop ####
 
 # create univ empty output table
 univ.table <- as.data.frame(matrix(ncol=8, nrow = (nrow(loci_analyse))))
 colnames(univ.table) <- c("LOC","CHR","START","STOP","phen","h2.obs","p","n_snps")
 univ.table$LOC = rep(loci_analyse$LOC)
 univ.table$CHR = rep(loci_analyse$CHR)
 univ.table$START = rep(loci_analyse$START)
 univ.table$STOP = rep(loci_analyse$STOP)
 
 i = 1 # initiate sub-loop ####
 while(i<nrow(loci_analyse)+1){ # nrow(loci_analyse)+1
  
  print(paste("Chromosome", chr, ":", "Processing locus",i,"of",nrow(loci_analyse), sep = " "))
  
  locus = process.locus(loci_analyse[i,], input)
  
  if(is.null(locus)==TRUE){
   i = i+1
  }else if(is.null(locus)==FALSE){
   
   univ.data = run.univ(locus)
   
   if(nrow(univ.data) == 0){
    i = i+1
   }
   
   if(nrow(univ.data) == 1){
    univ.table[univ.table$LOC == loci_analyse[i,]$LOC,]$phen[1] = univ.data$phen
    univ.table[univ.table$LOC == loci_analyse[i,]$LOC,]$h2.obs[1] = univ.data$h2.obs
    univ.table[univ.table$LOC == loci_analyse[i,]$LOC,]$p[1] = univ.data$p
    univ.table[univ.table$LOC == loci_analyse[i,]$LOC,]$n_snps[1] = locus$n.snps
    
    i = i+1
    
   }
   
  }
 }
 
 write.table(univ.table,
             file = paste(phenotype, "_", "univ_table_chr_", chr, ".txt", sep = ""),
             quote = FALSE,,
             row.names = FALSE)
 
 chr = chr + 1
 
}

# List all files in the working directory that contain "univ_table_chr_"
files_to_merge <- list.files(pattern = "univ_table_chr_")

univ_merged = fread(files_to_merge[1])

i = 2
while(i<length(files_to_merge)+1){
 
 print(paste("merging output for chromosome: ", i, sep = ""))
 
 add_this = fread(files_to_merge[i])
 univ_merged = rbind(univ_merged,
                     add_this)
 i = i+1
 
}

# Get the current timestamp
timestamp <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")

write.table(univ_merged,
            file = paste(phenotype,
                         "LAVA_univ",
                         timestamp,
                         ".txt"),
            quote = FALSE,
            row.names = FALSE)
