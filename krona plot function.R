library(dplyr)
library(tidyr)
library(readr)
library(phyloseq)

# -----------------------------
# 1️⃣ Ensure sample IDs exist
sample_data(ps1_agg)$dupl_id <- rownames(sample_data(ps1_agg))

# -----------------------------
# 2️⃣ Extract OTU table, taxonomy, metadata
otu <- as.data.frame(otu_table(ps1_agg))
tax <- as.data.frame(tax_table(ps1_agg))
meta <- as.data.frame(sample_data(ps1_agg))

# Transpose if OTUs are in columns
if(!taxa_are_rows(ps1_agg)){
  otu <- t(otu) %>% as.data.frame()
}

# Add OTU IDs
otu <- otu %>% rownames_to_column("OTU_ID")
tax <- tax %>%
  rownames_to_column("OTU_ID") %>%
  mutate(across(Kingdom:Species_exact, ~ifelse(is.na(.), "unclassified", .))) %>%
  unite("Taxonomy", Kingdom:Species_exact, sep = ";")

# -----------------------------
# 3️⃣ Merge OTU counts and taxonomy
otu_tax <- otu %>%
  left_join(tax, by = "OTU_ID")

# Identify sample columns
sample_cols <- intersect(colnames(otu_tax), rownames(meta))

# Pivot longer over samples
otu_long <- otu_tax %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "dupl_id",
    values_to = "Count"
  )

# Merge metadata
otu_meta <- otu_long %>%
  left_join(meta, by = "dupl_id")

# -----------------------------
# 4️⃣ Aggregate counts by category & taxonomy
grouped <- otu_meta %>%
  group_by(category, Taxonomy) %>%
  summarise(Total = sum(Count), .groups = "drop")

# Pivot so each group is a column
krona_input <- grouped %>%
  pivot_wider(names_from = category, values_from = Total, values_fill = 0)

# Write full group table (optional)
write_tsv(krona_input, "krona_input_by_group.txt")

# -----------------------------
# 5️⃣ Create per-category files + Krona HTML
krona_perl_path <- "C:/KronaTools/bin/ktImportText.pl"

unique_categories <- unique(grouped$category)
for(cat in unique_categories){
  
  cat_data <- otu_meta %>%
    filter(category == cat) %>%
    group_by(Taxonomy) %>%
    summarise(Total = sum(Count), .groups = "drop") %>%
    filter(Total > 0) %>%
    mutate(Taxonomy = gsub(";+",";", Taxonomy),
           Taxonomy = gsub(";+$","", Taxonomy))
  
  # Skip empty groups
  if(nrow(cat_data) == 0) next
  
  # Add "Other <1%" category
  total_counts <- sum(cat_data$Total)
  cat_data <- cat_data %>%
    mutate(RelAbund = Total / total_counts * 100)
  
  main_taxa <- cat_data %>%
    filter(RelAbund >= 1) %>%
    dplyr::select(Total, Taxonomy)
  
  other_taxa <- cat_data %>%
    filter(RelAbund < 1) %>%
    summarise(Total = sum(Total)) %>%
    mutate(Taxonomy = "Other") %>%
    dplyr::select(Total, Taxonomy)
  
  krona_cleaned <- bind_rows(main_taxa, other_taxa)
  
  # Write file (numeric first, taxonomy second)
  out_file <- paste0("krona_input_", cat, ".txt")
  write_tsv(krona_cleaned, out_file, col_names = FALSE)
  
  # Generate Krona HTML
  output_html <- paste0("krona_plot_", cat, ".html")
  cmd <- sprintf('perl "%s" -o "%s" "%s"', krona_perl_path, output_html, out_file)
  system(cmd)
  
  message("Krona plot generated: ", output_html)
}

# -----------------------------
# 6️⃣ Troubleshoot: inspect one file
read_tsv(out_file, col_names = FALSE)




























generate_krona <- function(ps_obj, category_col = "category", file_prefix = "krona") {
  library(dplyr)
  library(tidyr)
  library(readr)
  
  # Ensure dupl_id exists
  sample_data(ps_obj)$dupl_id <- rownames(sample_data(ps_obj))
  
  # OTU, taxonomy, metadata
  otu <- as.data.frame(otu_table(ps_obj))
  tax <- as.data.frame(tax_table(ps_obj))
  meta <- as.data.frame(sample_data(ps_obj))
  
  # Transpose OTU table if needed
  if(!taxa_are_rows(ps_obj)){
    otu <- t(otu)
    otu <- as.data.frame(otu)
  }
  
  # Add OTU IDs
  otu <- otu %>% rownames_to_column("OTU_ID")
  tax <- tax %>%
    rownames_to_column("OTU_ID") %>%
    mutate(across(Kingdom:Species_exact, ~ifelse(is.na(.), "unclassified", .))) %>%
    unite("Taxonomy", Kingdom:Species_exact, sep = ";")
  
  # Combine OTU + taxonomy
  otu_tax <- otu %>% left_join(tax, by = "OTU_ID")
  
  # Identify sample columns
  sample_cols <- intersect(colnames(otu_tax), rownames(meta))
  
  # Pivot longer
  otu_long <- otu_tax %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "dupl_id",
      values_to = "Count"
    )
  
  # Merge with metadata
  otu_meta <- otu_long %>% left_join(meta, by = "dupl_id")
  
  # Sum counts by group
  grouped <- otu_meta %>%
    group_by(.data[[category_col]], Taxonomy) %>%
    summarise(Total = sum(Count), .groups = "drop")
  
  # Write Krona-ready files per category
  unique_groups <- unique(grouped[[category_col]])
  krona_perl_path <- "C:/KronaTools/bin/ktImportText.pl"  # adjust path if needed
  
  for (grp in unique_groups) {
    cat_data <- otu_meta %>%
      filter(.data[[category_col]] == grp) %>%
      group_by(Taxonomy) %>%
      summarise(Total = sum(Count), .groups = "drop") %>%
      filter(Total > 0) %>%
      mutate(Taxonomy = gsub(";+",";", Taxonomy),
             Taxonomy = gsub(";+$","", Taxonomy))
    
    if(nrow(cat_data) == 0) next
    
    out_file <- paste0(file_prefix, "_", grp, ".txt")
    write_tsv(dplyr::select(cat_data, Total, Taxonomy), out_file, col_names = FALSE)
    
    output_html <- paste0(file_prefix, "_", grp, ".html")
    cmd <- sprintf('perl "%s" -o "%s" "%s"', krona_perl_path, output_html, out_file)
    system(cmd)
    message("Krona plot generated: ", output_html)
  }
}
generate_krona(ps1_agg, category_col = "category", file_prefix = "ps_filtered")
generate_krona(ps1s_agg, category_col = "category", file_prefix = "ps_ext_filtered")

####Troubleshooting#####
#Inspect file 
# Read the cleaned file without headers
krona_check <- read_tsv("krona_input_cleaned.txt", col_names = FALSE)
head(krona_check, 10)
str(krona_check)
