## ============================================================
## Drug Target Annotation Pipeline: Open Targets + OmniPath
## ============================================================

# Load libraries ####
library(otargen)
library(biomaRt)
library(data.table)
library(dplyr)
library(OmnipathR)

# Output directory ####
out_dir <- "/output/dir"
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

## ============================================================
## Step 1. Import OmniPath interaction network
## ============================================================
omnipath_interactions <- import_omnipath_interactions()
omnipath_genes <- unique(c(omnipath_interactions$source, omnipath_interactions$target))

## ============================================================
## Step 2. Import molecular exposures (UKB-PPC cis-pQTLs etc.)
## ============================================================
exposure_file <- "/path/to/LDclmpd_cispQTL_UKBBPPP.txt"
exposures <- fread(exposure_file)
protein_exposures <- unique(exposures$id.exposure)

## Combine exposures + interactors
instrument_universe <- unique(c(protein_exposures, omnipath_genes))

## ============================================================
## Step 3. Map to Ensembl Gene IDs (biomaRt)
## ============================================================
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(
 attributes = c("hgnc_symbol", "ensembl_gene_id"),
 filters = "hgnc_symbol",
 values = instrument_universe,
 mart = mart
)

saveRDS(mapping,
        file = file.path(out_dir, paste0("mapping_", timestamp, ".rda")))

## ============================================================
## Step 4. Query Open Targets for known drugs
## ============================================================
open_targets_annotation <- data.frame()

for (i in seq_len(nrow(mapping))) {
 message("Querying protein #", i, " of ", nrow(mapping))
 
 query_res <- knownDrugsQuery(mapping$ensembl_gene_id[i],
                              cursor = NULL,
                              freeTextQuery = NULL,
                              size = 10)
 
 if (!is.null(query_res)) {
  colnames(query_res) <- paste("OPENTARGETS", colnames(query_res), sep = "_")
  query_res$query_protein <- mapping$hgnc_symbol[i]
  query_res$query_protein_ensid <- mapping$ensembl_gene_id[i]
  
  open_targets_annotation <- rbind(open_targets_annotation, query_res)
 }
}

saveRDS(open_targets_annotation,
        file = file.path(out_dir, paste0("open_targets_annotation_", timestamp, ".rda")))

## ============================================================
## Step 5. Parse drug mechanisms (MOA)
## ============================================================
# Extract MOA type and targets
MOA_actionType <- sapply(open_targets_annotation$OPENTARGETS_drug.mechanismsOfAction.rows, function(x) {
 if (!is.null(x) && is.data.frame(x) && nrow(x) > 0) x$actionType[1] else NA
})

open_targets_annotation$MOA_actionType <- MOA_actionType

# Map to MR-interpretable categories
moa_map <- c(
 "INHIBITOR" = "inhibits", "ANTAGONIST" = "inhibits",
 "ANTISENSE INHIBITOR" = "inhibits", "RNAI INHIBITOR" = "inhibits",
 "BLOCKER" = "inhibits", "NEGATIVE MODULATOR" = "inhibits",
 "NEGATIVE ALLOSTERIC MODULATOR" = "inhibits", "INVERSE AGONIST" = "inhibits",
 
 "AGONIST" = "activates", "PARTIAL AGONIST" = "activates",
 "ACTIVATOR" = "activates", "OPENER" = "activates",
 "RELEASING AGENT" = "activates", "POSITIVE MODULATOR" = "activates",
 "POSITIVE ALLOSTERIC MODULATOR" = "activates", "STABILISER" = "activates",
 
 "EXOGENOUS PROTEIN" = "mimics", "EXOGENOUS GENE" = "mimics",
 "SUBSTRATE" = "mimics",
 
 "PROTEOLYTIC ENZYME" = "ambiguous", "HYDROLYTIC ENZYME" = "ambiguous",
 "DISRUPTING AGENT" = "ambiguous", "CROSS-LINKING AGENT" = "ambiguous",
 "BINDING AGENT" = "ambiguous", "MODULATOR" = "ambiguous", "OTHER" = "ambiguous",
 
 "VACCINE ANTIGEN" = "vaccine"
)

open_targets_annotation$effect_on_target <- moa_map[open_targets_annotation$MOA_actionType]

saveRDS(open_targets_annotation,
        file = file.path(out_dir, paste0("opentargets_annotated_", timestamp, ".rda")))

## ============================================================
## Step 6. Harmonise Open Targets indications to study phenotypes
## ============================================================

# Path to curated phenotype mapping CSV (EDIT PATH IF NEEDED)
phenomap_csv <- "path/to/annotated_opentargets_diseases_v0-1.csv"

# Load mapping and keep only rows flagged as relevant_phenotype == TRUE
pheno_map <- read.csv(phenomap_csv, check.names = FALSE)
pheno_map <- pheno_map[pheno_map$relevant_phenotype == TRUE, ]

# Join by Open Targets disease name
# (Assumes your Open Targets table has a column named "OPENTARGETS_disease.name")
open_targets_annotation <- open_targets_annotation %>%
 left_join(pheno_map, by = "OPENTARGETS_disease.name")

# Save the harmonised annotation
saveRDS(
 open_targets_annotation,
 file = file.path(out_dir, paste0("opentargets_annotated_harmonised_", timestamp, ".rda"))
)

