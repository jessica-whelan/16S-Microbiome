# Load required packages
library(tidyverse)  
library(dplyr)
library(tibble)
library(ggplot2)
library(readr)
library(phyloseq)
library(vegan)
library(ape)
library(ggpubr)
library(phangorn)
library(gghighlight)
library(pairwise)
library(DECIPHER)
library(mixOmics)
library(Biostrings)
library(pairwiseAdonis)
library(dada2)
library(gridExtra)
library(miLineage)
library(readxl)         # Excel import
library(decontam)
library(scales)
library(car)
library(openxlsx)
library(rgl)

# Import phyloseq object
physeq <- read_rds("dada2_phyloseq.rds")

# Extract OTU and taxonomy tables
otu <- otu_table(physeq)
tax <- tax_table(physeq)

# Helper function to load sample metadata
load_sample_data <- function(filepath) {
  read_excel(filepath) %>%
    column_to_rownames("sample") %>%
    sample_data()
}

# Import custom sample data
samples <- load_sample_data("sam_data excel for r.xlsx")
samples_noc <- load_sample_data("sam_data excel for r no c.xlsx")   #this is sample data without control mock community samples, totally optional

# Create phyloseq objects with different sample data
physeq <- phyloseq(otu, tax, samples)
physeq_noc <- phyloseq(otu, tax, samples_noc) #for sample data without mock community samples

# --- Decontam workflow (run before filtering) ---

# Define known mock taxa (we are inserting the genus names provided from Zymo mock community)
known_mock_taxa <- c(
  "Pseudomonas", "Escherichia-Shigella", "Salmonella",
  "Enterococcus", "Limosilactobacillus", "Staphylococcus",
  "Listeria", "Bacillus"
)

# Identify mock ASVs by genus
tax_df <- as.data.frame(tax_table(physeq)) %>%
  rownames_to_column("ASV_ID")

mock_asvs <- tax_df %>%
  filter(Genus %in% known_mock_taxa) %>%
  pull(ASV_ID)

# Identify control samples (category == "C")
control_samples <- sample_names(subset_samples(physeq, category == "C"))

# Access OTU table (ensure taxa are rows)
otu_mat <- otu_table(physeq)
if (!taxa_are_rows(physeq)) otu_mat <- t(otu_mat)

# Zero out mock ASVs in control samples
otu_mat[mock_asvs, control_samples] <- 0
otu_table(physeq) <- otu_mat

# Copy object for filtering
physeq_filtered <- physeq

# Flag negative controls
sample_data(physeq_filtered)$is_neg <-
  sample_data(physeq_filtered)$category == "C"

# Identify contaminants (prevalence method)
contam_df <- isContaminant(
  physeq_filtered,
  method = "prevalence",
  neg = "is_neg"
)

# Inspect results
table(contam_df$contaminant)
head(contam_df[contam_df$contaminant, ])

# Merge taxonomy for identified contaminants
tax_table_df <- as.data.frame(tax_table(physeq_filtered)) %>%
  rownames_to_column("ASV")

contam_df <- contam_df %>%
  rownames_to_column("ASV")

contam_with_tax <- left_join(contam_df, tax_table_df, by = "ASV")

# Extract contaminant ASVs
contaminant_asvs <- contam_with_tax %>%
  filter(contaminant) %>%
  pull(ASV)

# Prune contaminants
physeq_cleaned <- prune_taxa(
  !taxa_names(physeq_filtered) %in% contaminant_asvs,
  physeq_filtered
)

# To view identified and removed contaminants
contam_table <- contam_with_tax %>%
  dplyr::filter(contaminant) %>%
  dplyr::select(ASV, Kingdom, Phylum, Class, Order, Family, Genus)
View(contam_table)  # to view the table in separate R tab

# Remove control samples (category C)
physeq_decon <- subset_samples(physeq_cleaned, category != "C")

# Check number of taxa at each stage
ntaxa(physeq_filtered)   # before pruning
ntaxa(physeq_cleaned)    # after contaminant removal
ntaxa(physeq_decon)      # after removing controls

# --- Filtering and Pre-processing ---

# Initialize tracking table
filter_summary <- tibble(
  step = "raw (after decontam)",
  taxa = ntaxa(physeq_decon),
  samples = nsamples(physeq_decon),
  reads = sum(sample_sums(physeq_decon))
)

# 1. Remove mitochondrial taxa #Please note that for DADA2 the reference databased used is SILVA 138--
# human DNA, including mitochondrial and nuclear, is not in database, therefore there is no Family of Mitochondria 
# but this may be useful for other databases
physeq_nomit <- subset_taxa(physeq_decon, Family != "Mitochondria" | is.na(Family))
filter_summary <- add_row(filter_summary,
                          step = "remove mitochondria",
                          taxa = ntaxa(physeq_nomit),
                          samples = nsamples(physeq_nomit),
                          reads = sum(sample_sums(physeq_nomit))
)

# 2. Keep ASVs with confidence ≥ 0.95 (if confidence column exists)
tax_df <- as.data.frame(tax_table(physeq_nomit)) %>%
  rownames_to_column("ASV")

if ("confidence" %in% colnames(tax_df)) {
  high_conf_asvs <- tax_df %>%
    filter(confidence >= 0.95) %>%
    pull(ASV)
  
  physeq_nomit <- prune_taxa(high_conf_asvs, physeq_nomit)
  filter_summary <- add_row(filter_summary,
                            step = "confidence ≥ 0.95",
                            taxa = ntaxa(physeq_nomit),
                            samples = nsamples(physeq_nomit),
                            reads = sum(sample_sums(physeq_nomit))
  )
} else {
  message("No 'confidence' column found in taxonomy table — skipping confidence filter.")
}

# 3. Remove taxa with missing or empty Phylum
physeq_nomit <- subset_taxa(physeq_nomit, !is.na(Phylum) & Phylum != "")
filter_summary <- add_row(filter_summary,
                          step = "remove missing/empty phylum",
                          taxa = ntaxa(physeq_nomit),
                          samples = nsamples(physeq_nomit),
                          reads = sum(sample_sums(physeq_nomit))
)

# 4. Remove taxa with zero counts across all samples
physeq_nohdna <- prune_taxa(taxa_sums(physeq_nomit) > 0, physeq_nomit)
filter_summary <- add_row(filter_summary,
                          step = "remove zero-count taxa",
                          taxa = ntaxa(physeq_nohdna),
                          samples = nsamples(physeq_nohdna),
                          reads = sum(sample_sums(physeq_nohdna))
)

# Final filtered object
deconps <- physeq_nohdna



#--------------------Mock Community Analysis (Only needs to be done once)
#FIRST go back to Decontam workflow and do not remove group C. Do not zero out known microbes.  Just apply the known contaminant to the phyloseq object and move forward.
# 1. Subset phyloseq object to only mock community samples (category = "C")
ps_mock <- subset_samples(deconps, category == "C")

# 2. Aggregate to Genus (can also use "Family" if you prefer)
ps_genus <- tax_glom(ps_mock, taxrank = "Genus", NArm = FALSE)

# 3. Transform to relative abundance
ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

# 4. Melt to dataframe
df_mock <- psmelt(ps_rel)

# 5. Summarize across all mock samples
genus_table <- df_mock %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  mutate(Percentage = round(100 * TotalAbundance / sum(TotalAbundance), 2)) %>%
  arrange(desc(Percentage))
# Optional- view the table
view(genus_table)
# 6. (Optional) Collapse low-abundance genera into "Other"
genus_table <- genus_table %>%
  mutate(Genus = ifelse(Percentage < 1, "Other", Genus)) %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(TotalAbundance), .groups = "drop") %>%
  mutate(Percentage = round(100 * TotalAbundance / sum(TotalAbundance), 2)) %>%
  arrange(desc(Percentage))

# 7. Plot pie chart
ggplot(genus_table, aes(x = "", y = Percentage, fill = Genus)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Genus, ": ", Percentage, "%")),
            position = position_stack(vjust = 0.5), size = 3) +
  labs(title = "Genus-Level Relative Abundance in Mock Community") +
  theme_void()

#Before moving to secondary filtering, go back and re-enter the dropoff of group C, and redo analyses/filtering up to the next point


# --- Secondary Filtering ---

# --- Prevalence Analysis ---
otu <- otu_table(deconps)
tax <- tax_table(deconps)
# Compute prevalence per ASV
prev <- apply(otu, MARGIN = ifelse(taxa_are_rows(deconps), 1, 2),
              function(x) sum(x > 0))
prev_df <- data.frame(Prevalence = prev, tax)
# Summarize prevalence by Phylum
prev_by_phylum <- prev_df %>%
  as_tibble(rownames = "ASV") %>%
  group_by(Phylum) %>%
  summarise(
    mean_prevalence = mean(Prevalence),
    sum_prevalence = sum(Prevalence),
    asv_count = n(),
    .groups = "drop"
  )
# Summarize prevalence by Class
prev_by_class <- prev_df %>%
  as_tibble(rownames = "ASV") %>%
  group_by(Class) %>%
  summarise(
    mean_prevalence = mean(Prevalence),
    sum_prevalence = sum(Prevalence),
    asv_count = n(),
    .groups = "drop"
  )
view(prev_by_phylum)

# Define phyla to exclude, single-count phyla as presumable crossread errors
filter_phyla <- c("Deferribacterota", "Euryarchaeota", "Spirochaetota", "Sva0485") #substitute with your single-count phyla
ps1 <- subset_taxa(deconps, !Phylum %in% filter_phyla)
filter_summary <- add_row(filter_summary,
                          step = "secondary filtering (exclude selected phyla)",
                          taxa = ntaxa(ps1),
                          samples = nsamples(ps1),
                          reads = sum(sample_sums(ps1))
)
# --- Check up on Filtering Summary ---
filter_summary

# --- Conservative vs Extensive Filtering ---
# Genera to filter out based on Salter et al. 2014 (common contaminants, kit/reagent and lab environment contaminants)
filter_genera <- c(
  "Afipia", "Delftia", "Acinetobacter", "Escherichia-Shigella", 
  "Methylobacterium-Methylorubrum", "Cutibacterium", "Sphingobacterium", 
  "Burkholderia-Caballeronia-Paraburkholderia", "Aquabacterium", "Asticcacaulis", 
  "Aurantimonas", "Beijerinckia", "Bosea", "Bradyrhizobium", "Brevundimonas", 
  "Caulobacter", "Craurococcus", "Devosia", "Hoefleae", "Mesorhizobium", 
  "Methylobacterium", "Novosphingobium", "Ochrobactrum", "Paracoccus", 
  "Pedomicrobium", "Phyllobacterium", "Rhizobium", "Roseomonas", "Sphingobium", 
  "Sphingomonas", "Sphingopyxis", "Acidovorax", "Azoarcus", "Azospira", 
  "Burkholderia", "Comamonas", "Cupriavidus", "Curvibacter", "Delftiae", 
  "Duganella", "Herbaspirillum", "Janthinobacterium", "Kingella", "Leptothrix", 
  "Limnobacter", "Massilia", "Methylophilus", "Methyloversatilis", "Oxalobacter", 
  "Pelomonas", "Polaromonas", "Ralstonia", "Schlegelella", "Sulfuritalea", 
  "Undibacterium", "Variovorax", "Acinetobactera", "Enhydrobacter", "Enterobacter", 
  "Escherichia", "Nevskia", "Pasteurella", "Pseudomonas", "Pseudoxanthomonas", 
  "Psychrobacter", "Stenotrophomonas", "Xanthomonas", "Gp2", "Aeromicrobium", 
  "Actinomyces", "Arthrobacter", "Beutenbergia", "Brevibacterium", 
  "Corynebacterium", "Curtobacterium", "Dietzia", "Geodermatophilus", 
  "Janibacter", "Kocuria", "Microbacterium", "Micrococcus", "Microlunatus", 
  "Patulibacter", "Propionibacterium", "Rhodococcus", "Tsukamurella", 
  "Chryseobacterium", "Dyadobacter", "Flavobacterium", "Hydrotalea", "Niastella", 
  "Olivibacter", "Parabacteroides", "Pedobacter", "Prevotella", "Wautersiella", 
  "Abiotrophia", "Bacillus", "Brevibacillus", "Brochothrix", "Facklamia", 
  "Lactobacillus", "Paenibacillus", "Ruminococcus", "Staphylococcus", 
  "Streptococcus", "Veillonella", "Fusobacterium"
)
# Fungal contaminant taxa, although these are not in the silva 138 database
# but may be useful for other databases
filter_contaminant_taxa <- c(
  "Candidatus Annandia", "Candidatus Obscuribacter",
  "Candidatus Protochlamydia", "Candidatus Thiophysa"
)

# Conservative filtering object: ps1 (already created earlier)
# Extensive filtering: remove extra genera + contaminants
ps1s <- ps1 %>%
  subset_taxa(!Genus %in% filter_genera) %>%
  subset_taxa(!Genus %in% filter_contaminant_taxa)

# Track filtering summary
filter_summary <- filter_summary %>%
  add_row(
    step = "extensive filtering (remove suspect genera + contaminants)",
    taxa = ntaxa(ps1s),
    samples = nsamples(ps1s),
    reads = sum(sample_sums(ps1s))
  )

filter_summary

# --- Post-Filtering Prevalence Analysis Function ---
compute_prevalence <- function(physeq_obj) {
  otu <- otu_table(physeq_obj)
  tax <- tax_table(physeq_obj)
  
  prev <- apply(
    otu,
    MARGIN = ifelse(taxa_are_rows(physeq_obj), 1, 2),
    function(x) sum(x > 0)
  )
  
  prev_df <- data.frame(Prevalence = prev, tax) %>%
    as_tibble(rownames = "ASV")
  
  list(
    by_phylum = prev_df %>%
      group_by(Phylum) %>%
      summarise(
        mean_prevalence = mean(Prevalence),
        sum_prevalence = sum(Prevalence),
        asv_count = n(),
        .groups = "drop"
      ),
    by_class = prev_df %>%
      group_by(Class) %>%
      summarise(
        mean_prevalence = mean(Prevalence),
        sum_prevalence = sum(Prevalence),
        asv_count = n(),
        .groups = "drop"
      )
  )
}

# Run prevalence summaries for ps1 (conservative) and ps1s (extensive)
prev_ps1 <- compute_prevalence(ps1)
prev_ps1s <- compute_prevalence(ps1s)

# Inspect summaries
prev_ps1$by_phylum
prev_ps1$by_class

prev_ps1s$by_phylum
prev_ps1s$by_class


##################### Aggregating Data #################
#we have techincal replicates for our samples varying between 1 (no technical replicates) and 4
#original samples are identified by "dupl_id"
aggregate_multiple_phyloseq <- function(ps_list, id_col = "dupl_id") {
  
  aggregate_single <- function(ps) {
    # OTU table as matrix
    otu_mat <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu_mat <- t(otu_mat)
    
    # Sample metadata
    sample_df <- as.data.frame(sample_data(ps), stringsAsFactors = FALSE, check.names = FALSE)
    sample_df$SampleID <- rownames(sample_df)
    sample_df[[id_col]] <- as.character(sample_df[[id_col]])
    
    # Combine OTU table and metadata
    otu_df <- as.data.frame(otu_mat)
    otu_df$SampleID <- rownames(otu_df)
    combined_df <- dplyr::left_join(otu_df, sample_df, by = "SampleID")
    
    # Aggregate OTU counts by sum
    otu_cols <- setdiff(colnames(otu_df), c("SampleID", colnames(sample_df)))
    otu_with_id <- combined_df %>%
      dplyr::select(all_of(c(id_col, otu_cols)))
    agg_otu_df <- otu_with_id %>%
      dplyr::group_by(.data[[id_col]]) %>%
      dplyr::summarise(across(all_of(otu_cols), sum, na.rm = TRUE),
                       .groups = "drop")
    
    # Aggregate metadata: first row per group
    agg_meta_df <- sample_df %>%
      dplyr::group_by(.data[[id_col]]) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      tibble::column_to_rownames(id_col)
    
    # Prepare OTU matrix
    otu_only <- agg_otu_df %>% dplyr::select(-all_of(id_col)) %>% as.data.frame()
    rownames(otu_only) <- agg_otu_df[[id_col]]
    otu_mat_new <- as.matrix(otu_only)
    
    # QC table
    dupl_ids <- rownames(agg_meta_df)
    n_reps <- sapply(dupl_ids, function(x) sum(sample_df[[id_col]] == x))
    total_counts_before <- sapply(dupl_ids, function(x) sum(rowSums(otu_mat[sample_df[[id_col]] == x, , drop = FALSE])))
    total_counts_after <- rowSums(otu_mat_new)
    
    qc_df <- data.frame(
      dupl_id = dupl_ids,
      n_replicates = n_reps,
      total_counts_before = total_counts_before,
      total_counts_after = total_counts_after,
      zero_before = total_counts_before == 0,
      zero_after = total_counts_after == 0
    )
    
    # Build new phyloseq object
    ps_new <- phyloseq(
      otu_table(otu_mat_new, taxa_are_rows = FALSE),
      sample_data(agg_meta_df),
      tax_table(ps)
    )
    
    return(list(phyloseq = ps_new, QC = qc_df))
  }
  
  # Apply aggregation to all phyloseq objects in the list
  result_list <- lapply(ps_list, aggregate_single)
  
  # Name the list if not already named
  if (is.null(names(result_list))) names(result_list) <- paste0("ps", seq_along(result_list))
  
  return(result_list)
}
physeq_objects <- list(conservative = ps1, extensive = ps1s)
agg_results <- aggregate_multiple_phyloseq(physeq_objects, id_col = "dupl_id")

# Access aggregated phyloseq objects
ps1_agg  <- agg_results$conservative$phyloseq
ps1s_agg <- agg_results$extensive$phyloseq

# Access QC tables
View(agg_results$conservative$QC)
View(agg_results$extensive$QC) 



################### Raw Abundance Plots #######################
plot_abundance_by_tax <- function(physeq_obj, tax_rank = "Phylum") {
  # Ensure dupl_id exists
  sample_data(physeq_obj)$dupl_id <- rownames(sample_data(physeq_obj))
  
  # Generate the barplot
  p <- plot_bar(physeq_obj, x = "dupl_id", y = "Abundance", fill = tax_rank) +
    scale_y_continuous(labels = scales::comma_format()) +
    geom_bar(aes(color = !!as.name(tax_rank), fill = !!as.name(tax_rank)), 
             stat = "identity", position = "stack") +
    facet_wrap(~category, scales = "free_x") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   plot.title = element_text(size = 14, face = "bold"))
  
  return(p)
}

# Example usage for ps1_agg conservative filtering
p_ps1_agg_phylum <- plot_abundance_by_tax(ps1_agg, tax_rank = "Phylum")+
  ggtitle("ps1_agg: Phylum-level")
p_ps1_agg_family <- plot_abundance_by_tax(ps1_agg, tax_rank = "Family")+
  ggtitle("ps1_agg: Family-level")

# Example usage for ps1s_agg extensive filtering
p_ps1s_agg_phylum <- plot_abundance_by_tax(ps1s_agg, tax_rank = "Phylum")+
  ggtitle("ps1s_agg: Phylum-level")
p_ps1s_agg_family <- plot_abundance_by_tax(ps1s_agg, tax_rank = "Family")+
  ggtitle("ps1s_agg: Family-level")

# Optionally, arrange plots together
grid.arrange(p_ps1_agg_phylum, p_ps1_agg_family,
             p_ps1s_agg_phylum, p_ps1s_agg_family, ncol = 2)
grid.arrange(p_ps1_agg_phylum,
             p_ps1s_agg_phylum, ncol = 1)

# --- Extensive filtering (ps1_agg) ---
p_extensive_phylum <- plot_abundance_by_tax(ps1s_agg, tax_rank = "Phylum") + ggtitle("Extensive filtering - Phylum")
p_extensive_family <- plot_abundance_by_tax(ps1s_agg, tax_rank = "Family") + ggtitle("Extensive filtering - Family")

grid.arrange(p_extensive_phylum, p_extensive_family, ncol = 1)

# --- Conservative filtering (ps1s_agg) ---
p_conservative_phylum <- plot_abundance_by_tax(ps1_agg, tax_rank = "Phylum") + ggtitle("Conservative filtering - Phylum")
p_conservative_family <- plot_abundance_by_tax(ps1_agg, tax_rank = "Family") + ggtitle("Conservative filtering - Family")

grid.arrange(p_conservative_phylum, p_conservative_family, ncol = 1)


# ===============================
# Modular Phyloseq UNIFRAC Workflow (time consuming, optional. Only if using UF and WUF beta div)
# ===============================

# -------------------------------
# 1. Build tree for a given phyloseq object
# -------------------------------
build_tree_for_physeq <- function(physeq_obj, fasta_file) {
  # Import ASV sequences
  asv_seqs <- readDNAStringSet(fasta_file)
  asv_seqs_sub <- asv_seqs[taxa_names(physeq_obj)]
  
  # Reconstruct phyloseq object with refseq
  physeq_obj <- phyloseq(
    otu_table(otu_table(physeq_obj), taxa_are_rows = TRUE),
    tax_table(tax_table(physeq_obj)),
    sample_data(sample_data(physeq_obj)),
    refseq = asv_seqs_sub
  )
  
  # Align sequences
  alignment <- AlignSeqs(asv_seqs_sub, anchor = NA)
  
  # Build phylogenetic tree
  phydat <- phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- dist.ml(phydat)
  treeNJ <- NJ(dm)
  fit <- pml(treeNJ, data = phydat)
  fitGTR <- optim.pml(fit, model = "GTR",
                      optInv = TRUE, optGamma = TRUE,
                      rearrangement = "stochastic")
  tree <- fitGTR$tree
  
  # Root tree at midpoint
  tree_rooted <- phangorn::midpoint(tree)
  
  # Merge tree into phyloseq
  physeq_tree <- merge_phyloseq(physeq_obj, tree_rooted)
  
  # Add dupl_id for plotting
  sample_data(physeq_tree)$dupl_id <- rownames(sample_data(physeq_tree))
  
  return(physeq_tree)
}

# -------------------------------
# 2. Run UniFrac distances
# -------------------------------
compute_unifrac <- function(physeq_obj) {
  dist_unifrac <- phyloseq::distance(physeq_obj, method = "unifrac")
  dist_wunifrac <- phyloseq::distance(physeq_obj, method = "wunifrac")
  return(list(UF = dist_unifrac, WUF = dist_wunifrac))
}

# Function to remove zero-read samples
remove_zero_samples <- function(physeq_obj) {
  zero_samples <- sample_names(physeq_obj)[sample_sums(physeq_obj) == 0]
  if(length(zero_samples) > 0){
    message("Removing zero-read samples: ", paste(zero_samples, collapse = ", "))
    physeq_obj <- prune_samples(!sample_names(physeq_obj) %in% zero_samples, physeq_obj)
  }
  return(physeq_obj)
}

# Prune zero-read samples
ps_conservative_tree_clean <- remove_zero_samples(ps_conservative_tree)
ps_extensive_tree_clean   <- remove_zero_samples(ps_extensive_tree)

# Now compute distances safely
dist_bray_conservative <- phyloseq::distance(ps_conservative_tree_clean, method = "bray")
dist_jaccard_conservative <- phyloseq::distance(ps_conservative_tree_clean, method = "jaccard")

dist_bray_extensive <- phyloseq::distance(ps_extensive_tree_clean, method = "bray")
dist_jaccard_extensive <- phyloseq::distance(ps_extensive_tree_clean, method = "jaccard")

# -------------------------------
# 3. Apply to both filtering schemes
# -------------------------------
fasta_file <- "ASV_seqs1.fasta"

# Conservative (ps1_agg)
ps_conservative_tree <- build_tree_for_physeq(ps1_agg, fasta_file)
dist_conservative <- compute_unifrac(ps_conservative_tree)

# Extensive (ps1s_agg)
ps_extensive_tree <- build_tree_for_physeq(ps1s_agg, fasta_file)
dist_extensive <- compute_unifrac(ps_extensive_tree)

# -------------------------------
# 4. Save trees for reuse
# -------------------------------
saveRDS(ps_conservative_tree, file = "ps_conservative_tree.rds")
saveRDS(ps_extensive_tree, file = "ps_extensive_tree.rds")

# ================================
# Load previously saved phyloseq objects with trees
# ================================
ps_conservative_tree <- readRDS("ps_conservative_tree.rds")
ps_extensive_tree <- readRDS("ps_extensive_tree.rds")
# -------------------------------
# Plot phylogenetic tree
# -------------------------------

plot_tree(ps_conservative_tree, color = "category", label.tips = "Genus", ladderize = "left")
plot_tree(ps_extensive_tree, color = "category", label.tips = "Genus", ladderize = "left")

# -------------------------------
# Phylogenetic tree plots, agglomerated
# -------------------------------
# Conservative tree
ps_genus_cons <- tax_glom(ps_conservative_tree, taxrank = "Genus")
phyloseq::plot_tree(
  ps_genus_cons, 
  color = "category", 
  shape = "Phylum", 
  label.tips = "Genus"
)

# Extensive tree
ps_genus_ext<- tax_glom(ps_extensive_tree, taxrank = "Genus")
phyloseq::plot_tree(
  ps_genus_ext, 
  color = "category", 
  shape = "dupl_id", 
  label.tips = "Genus"
)



  
# -------------------------------Plot rel abundance======================
plot_rel_abundance <- function(physeq_obj, tax_rank = "Phylum", min_abund = NULL, title = NULL) {
  
  # Identify zero-read samples
  zero_samples <- sample_names(physeq_obj)[sample_sums(physeq_obj) == 0]
  
  if(length(zero_samples) > 0){
    message("Removing zero-read samples: ", paste(zero_samples, collapse = ", "))
    physeq_obj <- prune_samples(!sample_names(physeq_obj) %in% zero_samples, physeq_obj)
  }
  
  # Transform to relative abundance
  physeq_rel <- transform_sample_counts(physeq_obj, function(x) x / sum(x))
  
  # Filter taxa by minimum abundance if specified
  if (!is.null(min_abund)) {
    mean_abund <- taxa_sums(physeq_rel) / nsamples(physeq_rel)
    keep_taxa <- names(mean_abund[mean_abund > min_abund])
    physeq_rel <- prune_taxa(keep_taxa, physeq_rel)
  }
  
  # Add dupl_id for plotting
  sample_data(physeq_rel)$dupl_id <- rownames(sample_data(physeq_rel))
  
  # Generate plot
  p <- plot_bar(physeq_rel, x = "dupl_id", y = "Abundance", fill = tax_rank) +
    scale_y_continuous(labels = scales::percent_format()) +
    geom_bar(aes_string(color = tax_rank, fill = tax_rank), stat = "identity", position = "stack") +
    facet_wrap(~category, scales = "free_x") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"))
  
  if (!is.null(title)) p <- p + ggtitle(title)
  
  # Return a list with plot and removed samples
  return(list(plot = p, removed_samples = zero_samples))
}

# -------------------------------
# Relative abundance plots
# -------------------------------
# Conservative
p_con_phylum <- plot_rel_abundance(ps_conservative_tree, tax_rank = "Phylum", title = "Conservative - Phylum")
p_con_family <- plot_rel_abundance(ps_conservative_tree, tax_rank = "Family", title = "Conservative - Family")
p_con_genus <- plot_rel_abundance(ps_conservative_tree, tax_rank = "Genus", min_abund = 0.01, title = "Conservative - Genus (>1%)")

# Extensive
p_ext_phylum <- plot_rel_abundance(ps_extensive_tree, tax_rank = "Phylum", title = "Extensive - Phylum")
p_ext_family <- plot_rel_abundance(ps_extensive_tree, tax_rank = "Family", title = "Extensive - Family")
p_ext_genus <- plot_rel_abundance(ps_extensive_tree, tax_rank = "Genus", min_abund = 0.01, title = "Extensive - Genus (>1%)")

# Arrange plots
grid.arrange(p_con_phylum$plot, p_con_family$plot, p_con_genus$plot, ncol = 1)
grid.arrange(p_ext_phylum$plot, p_ext_family$plot, p_ext_genus$plot, ncol = 1)
grid.arrange(p_con_phylum$plot, p_con_genus$plot, ncol = 1)
grid.arrange(p_ext_phylum$plot, p_ext_genus$plot, ncol = 1)
#check which samples were removed
p_con_phylum$removed_samples
p_ext_phylum$removed_samples #sample 4 was removed 


#############Alpha Diversity#######################
#Inspect Libraries
# Conservative filtering
lib.size_cons <- data.frame(Lib.size = rowSums(otu_table(ps_conservative_tree)),
                            Sample = rownames(otu_table(ps_conservative_tree)))
lib.size_cons$Sample <- factor(lib.size_cons$Sample, levels = lib.size_cons$Sample[order(lib.size_cons$Lib.size)])

ggplot(lib.size_cons, aes(x = Sample, y = Lib.size)) +
  geom_bar(stat = "identity") +
  ggpubr::rotate_x_text() +
  xlab("Sample") + ylab("Library Size") +
  scale_y_continuous(labels = scales::comma_format()) +
  gghighlight(Lib.size < 1370, label_key = Sample) #based on your histogram highlight relevant lib size

# Extensive filtering
lib.size_ext <- data.frame(Lib.size = rowSums(otu_table(ps_extensive_tree)),
                           Sample = rownames(otu_table(ps_extensive_tree)))
lib.size_ext$Sample <- factor(lib.size_ext$Sample, levels = lib.size_ext$Sample[order(lib.size_ext$Lib.size)])

ggplot(lib.size_ext, aes(x = Sample, y = Lib.size)) +
  geom_bar(stat = "identity") +
  ggpubr::rotate_x_text() +
  xlab("Sample") + ylab("Library Size") +
  scale_y_continuous(labels = scales::comma_format()) +
  gghighlight(Lib.size < 823, label_key = Sample) #based on your histogram highlight relevant lib size

#prepare OTU tables
# Conservative
tab_cons <- as(otu_table(ps_conservative_tree), "matrix")
if(taxa_are_rows(ps_conservative_tree)) tab_cons <- t(tab_cons)
tab_cons <- tab_cons[rowSums(tab_cons) > 0, ]

# Extensive
tab_ext <- as(otu_table(ps_extensive_tree), "matrix")
if(taxa_are_rows(ps_extensive_tree)) tab_ext <- t(tab_ext)
tab_ext <- tab_ext[rowSums(tab_ext) > 0, ]

# Define a 2x2 layout
layout(matrix(1:4, ncol=2, byrow=TRUE))
#rarefaction curves
# Conservative
rarecurve(tab_cons, step=25, lwd=2, col="blue", ylab="Species Count", xlab="Library Size", label=F, cex=0.5)

# Extensive
rarecurve(tab_ext, step=25, lwd=2, col="red", ylab="Species Count", xlab="Library Size", label=F, cex=0.5)

#Rarefied vs observed species
# Conservative
raremax_cons <- min(rowSums(tab_cons))
Srare_cons <- rarefy(tab_cons, raremax_cons)
S_cons <- specnumber(tab_cons)
plot(S_cons, Srare_cons, xlab="Observed No. of Species", ylab="Rarefied No. of Species")
abline(0,1)

# Extensive
raremax_ext <- min(rowSums(tab_ext))
Srare_ext <- rarefy(tab_ext, raremax_ext)
S_ext <- specnumber(tab_ext)
plot(S_ext, Srare_ext, xlab="Observed No. of Species", ylab="Rarefied No. of Species")
abline(0,1)
#you should see a 2x2 plot now
#Rarefaction slope per sample
# Conservative
slopes_cons <- apply(tab_cons, 1, function(x) rareslope(x, sample=raremax_cons))
slopes_cons

# Extensive
slopes_ext <- apply(tab_ext, 1, function(x) rareslope(x, sample=raremax_ext))
slopes_ext

#Rarefy object
set.seed(708)

# Conservative (1370)
physeq.rarefied_cons <- rarefy_even_depth(ps_conservative_tree, sample.size = 1370, rngseed = 1, replace=FALSE) #change rarefied value as needed
min(sample_sums(physeq.rarefied_cons))

# Extensive (823)
physeq.rarefied_ext <- rarefy_even_depth(ps_extensive_tree, sample.size = 823, rngseed = 1, replace=FALSE) #change rarefied value as needed
min(sample_sums(physeq.rarefied_ext))


#Alpha Div Metrics

# Conservative
alpha_div_cons <- estimate_richness(physeq.rarefied_cons, measures=c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson"))
#optional: 
write.xlsx(alpha_div_cons, file="alpha_div_conservative.xlsx", rowNames=TRUE)

# Extensive
alpha_div_ext <- estimate_richness(physeq.rarefied_ext, measures=c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson"))
#optional:
write.xlsx(alpha_div_ext, file="alpha_div_extensive.xlsx", rowNames=TRUE)

# Combine metadata
alpha_div.df_cons <- cbind(as.data.frame(sample_data(physeq.rarefied_cons))[,1:3], alpha_div_cons)
alpha_div.df_ext <- cbind(as.data.frame(sample_data(physeq.rarefied_ext))[,1:3], alpha_div_ext)

# Normality
shapiro.test(alpha_div_cons$Shannon)
shapiro.test(alpha_div_ext$Shannon)
# T-tests or Wilcoxon (adjust if non-normal)


#Plot Alpha Div
# Conservative plots
# Chao1
P_chao1_cons <- plot_richness(physeq.rarefied_cons, x="category", measures="Chao1", color="category") +
  geom_boxplot(alpha=0.6, outlier.colour="red", outlier.shape=8)  +
  stat_compare_means(method="t.test", label="p.format") +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))

# Shannon
P_shannon_cons <- plot_richness(physeq.rarefied_cons, x="category", measures="Shannon", color="category") +
  geom_boxplot(alpha=0.6, outlier.colour="red", outlier.shape=8) +  
  stat_compare_means(method="t.test", label="p.format") +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))

# InvSimpson
P_invsimpson_cons <- plot_richness(physeq.rarefied_cons, x="category", measures="InvSimpson", color="category") +
  geom_boxplot(alpha=0.6, outlier.colour="red", outlier.shape=8) +  
  stat_compare_means(method="t.test", label="p.format") +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))

# Combine into a single row (Conservative only)
ggarrange(P_chao1_cons, P_shannon_cons, P_invsimpson_cons,
          ncol=3, nrow=1)
          
#Plot Extensive Filtering
# Chao1
P_chao1_ext <- plot_richness(physeq.rarefied_ext, x="category", measures="Chao1", color="category") +
  geom_boxplot(alpha=0.6, outlier.colour="red", outlier.shape=8) +  
  stat_compare_means(method="wilcox.test", label="p.format") +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))

# Shannon
P_shannon_ext <- plot_richness(physeq.rarefied_ext, x="category", measures="Shannon", color="category") +
  geom_boxplot(alpha=0.6, outlier.colour="red", outlier.shape=8) +  
  stat_compare_means(method="wilcox.test", label="p.format") +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))

# InvSimpson
P_invsimpson_ext <- plot_richness(physeq.rarefied_ext, x="category", measures="InvSimpson", color="category") +
  geom_boxplot(alpha=0.6, outlier.colour="red", outlier.shape=8) +  
  stat_compare_means(method="wilcox.test", label="p.format") +
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))

# Combine into a single row (Extensive only)
ggarrange(P_chao1_ext, P_shannon_ext, P_invsimpson_ext,
          ncol=3, nrow=1 )
        

# Combine all plots into a single figure, not optimized, titles overlap
ggarrange(
  P_chao1_cons, P_shannon_cons, P_invsimpson_cons,
  P_chao1_ext, P_shannon_ext, P_invsimpson_ext,
  ncol=3, nrow=2,
  labels=c("Chao1 - Conservative", "Shannon - Conservative", "InvSimpson - Conservative",
           "Chao1 - Extensive", "Shannon - Extensive", "InvSimpson - Extensive")
)


######################### Modular Beta Diversity Workflow #########################
# -------------------------------
# Function to run beta diversity for multiple metrics
# -------------------------------
#########################
# Beta Diversity Workflow
#########################
beta_diversity_analysis <- function(physeq_obj, 
                                    category_col = "category", 
                                    title_prefix = "Beta Diversity",
                                    metrics = c("bray", "jaccard", "unifrac", "wunifrac")) {
  
  # Remove zero-read samples
  zero_samples <- sample_names(physeq_obj)[sample_sums(physeq_obj) == 0]
  if(length(zero_samples) > 0){
    message("Removing zero-read samples: ", paste(zero_samples, collapse = ", "))
    physeq_obj <- prune_samples(!sample_names(physeq_obj) %in% zero_samples, physeq_obj)
  }
  
  # Normalize counts to relative abundance (for Bray & Jaccard)
  physeq_rel <- transform_sample_counts(physeq_obj, function(x) {
    if(sum(x) > 0) x / sum(x) else x
  })
  
  # Extract metadata
  metadata <- as(sample_data(physeq_obj), "data.frame")
  
  results <- list()
  
  # ---- Bray-Curtis ----
  if("bray" %in% metrics){
    distBC <- phyloseq::distance(physeq_rel, method = "bray")
    set.seed(42)
    adonisBC <- adonis2(distBC ~ metadata[[category_col]])
    betadisperBC <- betadisper(distBC, metadata[[category_col]])
    ordBC <- ordinate(physeq_rel, method = "PCoA", distance = distBC)
    pBC <- plot_ordination(physeq_rel, ordBC, color = category_col, shape = category_col) +
      geom_point(size = 3) +
      stat_ellipse() +
      annotate("text", x = 0, y = 0.3, 
               label = paste("PERMANOVA p =", formatC(adonisBC$`Pr(>F)`[1], format = "f", digits = 3))) +
      ggtitle(paste(title_prefix, "- Bray-Curtis PCoA")) +
      theme(text = element_text(size = 16))
    
    results$bray <- list(distance = distBC, permanova = adonisBC, betadisper = betadisperBC, plot = pBC)
  }
  
  # ---- Jaccard ----
  if("jaccard" %in% metrics){
    distJ <- phyloseq::distance(physeq_rel, method = "jaccard")
    set.seed(42)
    adonisJ <- adonis2(distJ ~ metadata[[category_col]])
    betadisperJ <- betadisper(distJ, metadata[[category_col]])
    ordJ <- ordinate(physeq_rel, method = "PCoA", distance = distJ)
    pJ <- plot_ordination(physeq_rel, ordJ, color = category_col, shape = category_col) +
      geom_point(size = 3) +
      stat_ellipse() +
      annotate("text", x = 0, y = 0.3, 
               label = paste("PERMANOVA p =", formatC(adonisJ$`Pr(>F)`[1], format = "f", digits = 3))) +
      ggtitle(paste(title_prefix, "- Jaccard PCoA")) +
      theme(text = element_text(size = 16))
    
    results$jaccard <- list(distance = distJ, permanova = adonisJ, betadisper = betadisperJ, plot = pJ)
  }
  
  # ---- UniFrac (requires tree) ----
  if("unifrac" %in% metrics){
    if(is.null(phy_tree(physeq_obj))){
      warning("UniFrac requested but phylogenetic tree not found. Skipping.")
    } else {
      distUF <- phyloseq::distance(physeq_obj, method = "unifrac")
      set.seed(42)
      adonisUF <- adonis2(distUF ~ metadata[[category_col]])
      betadisperUF <- betadisper(distUF, metadata[[category_col]])
      ordUF <- ordinate(physeq_obj, method = "PCoA", distance = distUF)
      pUF <- plot_ordination(physeq_obj, ordUF, color = category_col, shape = category_col) +
        geom_point(size = 3) +
        stat_ellipse() +
        annotate("text", x = 0, y = 0.3, 
                 label = paste("PERMANOVA p =", formatC(adonisUF$`Pr(>F)`[1], format = "f", digits = 3))) +
        ggtitle(paste(title_prefix, "- UniFrac PCoA")) +
        theme(text = element_text(size = 16))
      
      results$unifrac <- list(distance = distUF, permanova = adonisUF, betadisper = betadisperUF, plot = pUF)
    }
  }
  
  # ---- Weighted UniFrac ----
  if("wunifrac" %in% metrics){
    if(is.null(phy_tree(physeq_obj))){
      warning("Weighted UniFrac requested but phylogenetic tree not found. Skipping.")
    } else {
      distWUF <- phyloseq::distance(physeq_obj, method = "wunifrac")
      set.seed(42)
      adonisWUF <- adonis2(distWUF ~ metadata[[category_col]])
      betadisperWUF <- betadisper(distWUF, metadata[[category_col]])
      ordWUF <- ordinate(physeq_obj, method = "PCoA", distance = distWUF)
      pWUF <- plot_ordination(physeq_obj, ordWUF, color = category_col, shape = category_col) +
        geom_point(size = 3) +
        stat_ellipse() +
        annotate("text", x = 0, y = 0.3, 
                 label = paste("PERMANOVA p =", formatC(adonisWUF$`Pr(>F)`[1], format = "f", digits = 3))) +
        ggtitle(paste(title_prefix, "- Weighted UniFrac PCoA")) +
        theme(text = element_text(size = 16))
      
      results$wunifrac <- list(distance = distWUF, permanova = adonisWUF, betadisper = betadisperWUF, plot = pWUF)
    }
  }
  
  return(results)
}

# Example usage

# Conservative
beta_cons <- beta_diversity_analysis(physeq.rarefied_cons, category_col = "category", title_prefix = "Conservative Filtering",
                                     metrics = c("bray","jaccard","unifrac","wunifrac"))
# Access plots
beta_cons$bray$plot
beta_cons$jaccard$plot
beta_cons$unifrac$plot
beta_cons$wunifrac$plot
#optionally plot together:
grid.arrange(beta_cons$bray$plot,
             beta_cons$jaccard$plot,
             beta_cons$unifrac$plot,
             beta_cons$wunifrac$plot,
             ncol=2)

# Print PERMANOVA p-values
cat("Conservative Filtering PERMANOVA p-values:\n")
cat("Bray-Curtis:", formatC(beta_cons$bray$permanova$`Pr(>F)`[1], digits = 3), "\n")
cat("Jaccard:", formatC(beta_cons$jaccard$permanova$`Pr(>F)`[1], digits = 3), "\n")
cat("UniFrac:", formatC(beta_cons$unifrac$permanova$`Pr(>F)`[1], digits = 3), "\n")
cat("Weighted UniFrac:", formatC(beta_cons$wunifrac$permanova$`Pr(>F)`[1], digits = 3), "\n\n")

# Extensive
beta_ext <- beta_diversity_analysis(physeq.rarefied_ext, category_col = "category", title_prefix = "Extensive Filtering",
                                    metrics = c("bray","jaccard","unifrac","wunifrac"))
# Access plots
beta_ext$bray$plot
beta_ext$jaccard$plot
beta_ext$unifrac$plot
beta_ext$wunifrac$plot
#optionally plot together:
grid.arrange(beta_ext$bray$plot,
             beta_ext$jaccard$plot,
             beta_ext$unifrac$plot,
             beta_ext$wunifrac$plot,
             ncol=2)
# Print PERMANOVA p-values
cat("Extensive Filtering PERMANOVA p-values:\n")
cat("Bray-Curtis:", formatC(beta_ext$bray$permanova$`Pr(>F)`[1], digits = 3), "\n")
cat("Jaccard:", formatC(beta_ext$jaccard$permanova$`Pr(>F)`[1], digits = 3), "\n")
cat("UniFrac:", formatC(beta_ext$unifrac$permanova$`Pr(>F)`[1], digits = 3), "\n")
cat("Weighted UniFrac:", formatC(beta_ext$wunifrac$permanova$`Pr(>F)`[1], digits = 3), "\n")



#########################SPLSDA###############################
# 1. Filtering with checks

library(rgl)  # for 3D plotting if needed

### ===================== ASV Filtering =====================
filter_asvs <- function(ps, prevalence = 0.05, abundance = 50, tax_rank = "Genus", plot = TRUE){
  prevdf <- apply(otu_table(ps), 
                  MARGIN = if(taxa_are_rows(ps)) 1 else 2,
                  FUN = function(x) sum(x > 0))
  
  prevdf <- data.frame(
    Prevalence = prevdf,
    TotalAbundance = taxa_sums(ps),
    tax_table(ps)
  )
  
  prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))
  
  if(plot){
    ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps), color=Phylum)) +
      geom_hline(yintercept = prevalence, linetype=2, alpha=0.5) +
      geom_point(size=2, alpha=0.7) +
      scale_x_log10() + xlab("Total Abundance") + ylab("Prevalence [Fraction Samples]") +
      facet_wrap(~Phylum) + theme(legend.position="none")
  }
  
  keepTaxa <- rownames(prevdf1)[
    (prevdf1$Prevalence >= prevalence * nsamples(ps)) |  #change between & and | as needed
      (prevdf1$TotalAbundance > abundance)
  ]
  
  ps_filtered <- prune_taxa(keepTaxa, ps)
  ps_agg <- tax_glom(ps_filtered, tax_rank, NArm=FALSE)
  
  cat("Original taxa:", ntaxa(ps), "\n")
  cat("After filtering:", ntaxa(ps_filtered), "\n")
  cat("After agglomeration:", ntaxa(ps_agg), "\n")
  
  return(ps_agg)
}

### ===================== CLR Preparation =====================
prepare_clr_data <- function(ps){
  OTU <- otu_table(ps)
  data.raw <- if(taxa_are_rows(ps)) as.data.frame(t(OTU)) else as.data.frame(OTU)
  
  # Remove zero-sum samples
  data.raw <- data.raw[rowSums(data.raw) > 0, ]
  
  # Normalize to relative abundance (TSS)
  data.rel <- sweep(data.raw, 1, rowSums(data.raw), FUN="/")
  
  # Add pseudocount
  data.pos <- data.rel + 1e-6
  
  # CLR transformation
  data.clr <- log(data.pos) - rowMeans(log(data.pos))
  
  # Sample metadata and taxonomy
  sample.meta <- as.data.frame(sample_data(ps))
  sample.meta <- sample.meta[rownames(data.clr), , drop=FALSE]
  taxonomy <- as.data.frame(tax_table(ps))
  
  # Checks
  stopifnot(!any(is.na(data.clr)))
  stopifnot(!any(data.clr == -Inf | data.clr == Inf))
  stopifnot(all(rownames(data.clr) == rownames(sample.meta)))
  
  list(data.raw=data.raw, data.clr=data.clr, taxonomy=taxonomy, sample.meta=sample.meta)
}

### ===================== PCA Plotting =====================
plot_pca <- function(data.clr, sample.meta, group_var="category", ind_names=TRUE){
  pca.res <- pca(data.clr, logratio="none")
  plotIndiv(pca.res, group=sample.meta[[group_var]], ind.names=ind_names,
            title="PCA Plot (CLR-transformed)", ellipse=TRUE, legend=TRUE)
  return(pca.res)
}

### ===================== Explore Optimal ncomp =====================
explore_ncomp <- function(data.clr, sample.meta, group_var="category", max_comp=5){
  Y <- factor(sample.meta[[group_var]])
  
  set.seed(123)
  perf.res <- perf(
    splsda(X = data.clr, Y = Y, ncomp = max_comp, keepX = rep(10, max_comp)),
    validation = "Mfold", folds = 5, nrepeat = 50, 
    dist = "centroids.dist", progressBar = TRUE
  )
  
  plot(perf.res, overlay="measure", sd=TRUE)
  perf.res
}

### ===================== Tune keepX =====================
tune_keepX <- function(data.clr, sample.meta, ncomp=2, group_var="category"){
  Y <- factor(sample.meta[[group_var]])
  
  grid.keepX <- seq(5, 100, 5)
  set.seed(123)
  tune.res <- tune.splsda(
    X = data.clr, Y = Y, ncomp = ncomp, logratio="none",
    test.keepX = grid.keepX, validation="Mfold", folds=5, nrepeat=50,
    dist="centroids.dist", measure="BER", progressBar=TRUE
  )
  
  optimal.keepX <- tune.res$choice.keepX[1:ncomp]
  list(tune.res=tune.res, optimal.keepX=optimal.keepX)
}

### ===================== Run sPLS-DA and Evaluate Performance =====================
run_splsda_with_perf <- function(ps, data.clr, sample.meta, ncomp=2, optimal.keepX=NULL, group_var="category"){
  Y <- factor(sample.meta[[group_var]])
  
  if(is.null(optimal.keepX)){
    keepX <- rep(10, ncomp)
  } else {
    keepX <- optimal.keepX
  }
  
  splsda.res <- splsda(data.clr, Y, ncomp=ncomp, keepX=keepX)
  
  # Plot samples
  plotIndiv(splsda.res, group=Y, ind.names=TRUE, ellipse=TRUE, legend=TRUE,
            title=paste("sPLS-DA: Comp 1 & 2"), style="ggplot2")
  
  # Loadings
  tax_table_df <- as.data.frame(tax_table(ps))
  tax_columns <- c("Class","Order","Family","Genus")
  short_tax_labels <- apply(tax_table_df[,tax_columns],1,function(x) paste(na.omit(x), collapse="; "))
  names(short_tax_labels) <- rownames(tax_table_df)
  
  for(comp in 1:ncomp){
    plotLoadings(splsda.res, comp=comp, method="median", contrib="max", ndisplay=10,
                 name.var=short_tax_labels, title=paste0("Loadings - Comp ", comp))
  }
  
  # Performance
  perf.res <- perf(splsda.res, validation="Mfold", folds=5, nrepeat=100, auc=TRUE, progressBar=TRUE)
  
  list(splsda.res=splsda.res, perf.res=perf.res)
}

### ===================== Optional 3D Plot =====================
plot_splsda_3d <- function(splsda.res, sample.meta, group_var="category"){
  Y <- factor(sample.meta[[group_var]])
  plotIndiv(splsda.res, comp=c(1,2,3), group=Y, ind.names=TRUE, style="3d", legend=TRUE,
            title="sPLSDA: 3D Plot of Components")
}
# ------------------ Conservative Filtering ------------------
ps <- ps1_agg   # or ps1s_agg for extensive filtering

# ------------------ 1. Filter ASVs ------------------
ps_filtered <- filter_asvs(ps, prevalence=0.05, abundance=50)

# ------------------ 2. Prepare CLR data ------------------
clr_data <- prepare_clr_data(ps_filtered)

# ------------------ 3. PCA ------------------
pca_res <- plot_pca(clr_data$data.clr, clr_data$sample.meta)

# ------------------ 4. Explore optimal ncomp ------------------
perf_ncomp <- explore_ncomp(clr_data$data.clr, clr_data$sample.meta, max_comp=5)

# ------------------ 5. Tune keepX ------------------
keepX_res <- tune_keepX(clr_data$data.clr, clr_data$sample.meta, ncomp=2) #enter optimal ncomp value

# ------------------ 6. Run sPLS-DA with optimized keepX and manual ncomp ------------------
splsda_final <- run_splsda_with_perf(ps_filtered, clr_data$data.clr, clr_data$sample.meta,
                                     ncomp=2, optimal.keepX=keepX_res$optimal.keepX) #enter optimal ncomp value

# Inspect error rates
splsda_final$perf.res$error.rate
splsda_final$perf.res$error.rate.class

# ------------------ 7. Optional 3D Plot (only for ncomp = 3 or more) ------------------
plot_splsda_3d(splsda_final$splsda.res, clr_data$sample.meta)

# Access optimized keepX
keepX_res$optimal.keepX




# ------------------ Extensive Filtering ------------------
ps_ext <- ps1s_agg   # Extensive filtering dataset

# 1. Filter ASVs
ps_ext_filtered <- filter_asvs(ps_ext, prevalence=0.05, abundance=30)  # Example thresholds for extensive filtering

# 2. Prepare CLR data
clr_data_ext <- prepare_clr_data(ps_ext_filtered)

# 3. PCA
pca_res_ext <- plot_pca(clr_data_ext$data.clr, clr_data_ext$sample.meta)

# 4. Explore optimal ncomp
perf_ncomp_ext <- explore_ncomp(clr_data_ext$data.clr, clr_data_ext$sample.meta, max_comp=5)

# 5. Tune keepX
keepX_res_ext <- tune_keepX(clr_data_ext$data.clr, clr_data_ext$sample.meta, ncomp=3) #enter optimal ncomp value

# 6. Run sPLS-DA with optimized keepX
splsda_final_ext <- run_splsda_with_perf(ps_ext_filtered, clr_data_ext$data.clr, clr_data_ext$sample.meta,
                                         ncomp=3, optimal.keepX=keepX_res_ext$optimal.keepX) #enter optimal ncomp value

# Inspect error rates
splsda_final_ext$perf.res$error.rate
splsda_final_ext$perf.res$error.rate.class

# 7. Optional 3D Plot
plot_splsda_3d(splsda_final_ext$splsda.res, clr_data_ext$sample.meta)

# Access optimized keepX
keepX_res_ext$optimal.keepX


