library(data.table)

protein_info_path = "path/to/9606.protein.info.v12.0.txt"

protein_links_path = "path/to/9606.protein.links.v12.0.txt"


protein_info_df = fread(protein_info_path)
protein_links_df = fread(protein_links_path)

str(protein_info_df)
str(protein_links_df)


# Add the protein1_name column by matching protein1 with string_protein_id in protein_info_df
protein_links_df[, protein1_name := protein_info_df[match(protein1, `#string_protein_id`), preferred_name]]

# Add the protein2_name column by matching protein2 with string_protein_id in protein_info_df
protein_links_df[, protein2_name := protein_info_df[match(protein2, `#string_protein_id`), preferred_name]]

# Generate the timestamp in the format YYYY_MM_DD_HH_MM_SS
timestamp <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")

# Create the filename with the timestamp
output_filename <- paste0("protein_links_annotated_genenames_", timestamp, ".txt")
output_filename = paste("output/directory",
                        output_filename,
                        sep = "")

# Save the updated protein_links_df to a .txt file with no strange delimiters
fwrite(protein_links_df, file = output_filename, sep = "\t", quote = FALSE, row.names = FALSE)