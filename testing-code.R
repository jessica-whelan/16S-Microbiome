

####testing new visuals###### STILL IN PROGRESS
#Heatmap on conservative filtering
short_tax_labels <- apply(tax_table(ps_filtered)[,c("Class","Order","Family","Genus")], 1, function(x) paste(na.omit(x), collapse="; "))
cim(t(splsda_final$splsda.res$X), comp = 1:2, row.names = short_tax_labels, margin = c(5,10))
#heatmap on extensive filtering
short_tax_labels <- apply(tax_table(ps_ext_filtered)[,c("Class","Order","Family","Genus")], 1, function(x) paste(na.omit(x), collapse="; "))
cim(t(splsda_final_ext$splsda.res$X), comp = 1:2, row.names = short_tax_labels, margin = c(5,10))


#Abundance plots
#-------------------------
# Helper functions
#-------------------------

# Prune taxa with zero counts
prune_zeros <- function(ps) {
  prune_taxa(taxa_sums(ps) > 0, ps)
}

# Aggregate phyloseq object to a taxonomic rank
aggregate_taxa <- function(ps, rank) {
  tax_glom(ps, taxrank = rank, NArm = FALSE)
}

# Transform counts to relative abundance
relative_abundance <- function(ps) {
  transform_sample_counts(ps, function(x) x / sum(x))
}

# Melt phyloseq object to dataframe, keep only non-zero abundances
ps_to_df <- function(ps) {
  df <- psmelt(ps)
  df <- df %>% filter(Abundance > 0)
  return(df)
}

# Filter top N taxa by abundance
top_n_taxa <- function(ps, N = 20) {
  top_taxa <- names(sort(taxa_sums(ps), decreasing = TRUE))[1:N]
  prune_taxa(top_taxa, ps)
}

#-------------------------
# Barplot function
#-------------------------
plot_bar_abundance <- function(ps, x = "dupl_id", fill = "Family", facet = "category", 
                               relative = FALSE, topN = NULL, label = FALSE) {
  
  ps <- prune_zeros(ps)
  
  # Ensure dupl_id exists
  if (!"dupl_id" %in% sample_variables(ps)) {
    sample_data(ps)$dupl_id <- rownames(sample_data(ps))
  }
  
  if (!is.null(topN)) {
    ps <- top_n_taxa(ps, topN)
  }
  
  if (relative) {
    ps <- relative_abundance(ps)
  }
  
  df <- ps_to_df(ps)
  
  # Check facet exists
  if (!(facet %in% colnames(df))) stop(paste("Facet column", facet, "not found in the phyloseq object"))
  
  p <- ggplot(df, aes_string(x = x, y = "Abundance", fill = fill)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +
    facet_wrap(as.formula(paste0("~", facet)), scales = "free_x") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  if (relative) p <- p + scale_y_continuous(labels = scales::percent_format())
  if (!relative) p <- p + scale_y_continuous(labels = scales::comma_format())
  
  if (label) {
    p <- p + geom_text(aes_string(label = fill), position = position_stack(vjust = 0.5), 
                       size = 2.5, check_overlap = TRUE)
  }
  
  return(p)
}

#-------------------------
# Piechart function
#-------------------------
plot_pie_abundance <- function(ps, rank = "Genus", facet = "category", topN = NULL, min_percent = 1) {
  
  ps <- prune_zeros(ps)
  ps <- aggregate_taxa(ps, rank)
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  ps <- relative_abundance(ps)
  
  # Ensure dupl_id exists (optional for labeling/facet)
  if (!"dupl_id" %in% sample_variables(ps)) {
    sample_data(ps)$dupl_id <- rownames(sample_data(ps))
  }
  
  df <- ps_to_df(ps)
  
  # Summarize relative abundance by facet and rank
  df_group <- df %>%
    group_by(across(all_of(c(facet, rank)))) %>%
    summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
    group_by(across(all_of(facet))) %>%
    mutate(RelAbundance = TotalAbundance / sum(TotalAbundance))
  
  # Collapse rare taxa (< min_percent) into "Other"
  df_group <- df_group %>%
    mutate(!!rank := ifelse(RelAbundance*100 < min_percent, "Other", !!sym(rank))) %>%
    group_by(across(all_of(c(facet, rank)))) %>%
    summarise(RelAbundance = sum(RelAbundance), .groups = "drop")
  
  # Keep top N taxa per facet if requested
  if (!is.null(topN)) {
    df_group <- df_group %>%
      group_by(across(all_of(facet))) %>%
      slice_max(order_by = RelAbundance, n = topN)
  }
  
  # Plot
  p <- ggplot(df_group, aes_string(x = '""', y = 'RelAbundance', fill = rank)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar("y") +
    facet_wrap(as.formula(paste0("~", facet))) +
    theme_void() +
    theme(legend.position = "right")
  
  return(p)
}

#-------------------------
# Example usage
#-------------------------
#conservative filtering
# Raw abundance barplot
p_raw <- plot_bar_abundance(ps_filtered, x = "dupl_id", fill = "Family", 
                            relative = FALSE, 
                            topN=40, #change this to whatever you want
                            label = TRUE)
p_raw
# Relative abundance barplot
p_rel <- plot_bar_abundance(ps_filtered, x = "dupl_id", topN=40,fill = "Family", relative = TRUE, label = TRUE)
p_rel
# Top 20 genus barplot
p_top20 <- plot_bar_abundance(ps_filtered, x = "dupl_id", fill = "Genus", relative = TRUE, topN = 40, label = TRUE)
p_top20
# Piechart for all samples
pie_all <- plot_pie_abundance(ps_filtered, rank = "Genus")
pie_all
# Piechart for mock community only not working currently
#pie_mock <- plot_pie_abundance(subset_samples(ps_filtered, category == "C"), rank = "Genus", topN = 20)

# Combine barplots
grid.arrange(nrow = 2, p_raw, p_rel)

####Extensive Filtering
# Raw abundance barplot
p_raws <- plot_bar_abundance(ps_ext_filtered, x = "dupl_id", fill = "Family", 
                             relative = FALSE, 
                             topN=40,
                             label = TRUE)
p_raws
# Relative abundance barplot
p_rels <- plot_bar_abundance(ps_ext_filtered, x = "dupl_id", topN=40,fill = "Family", relative = TRUE, label = TRUE)
p_rels
# Top 20 genus barplot
p_top20s <- plot_bar_abundance(ps_ext_filtered, x = "dupl_id",topN=40, fill = "Genus", relative = TRUE, label = TRUE)
p_top20s
# Piechart for all samples
pie_alls <- plot_pie_abundance(ps_ext_filtered, rank = "Family")
pie_alls
# Piechart for mock community only not working currently
#pie_mock <- plot_pie_abundance(subset_samples(ps1_agg, category == "C"), rank = "Genus", topN = 20)

# Combine barplots
grid.arrange(nrow = 2, p_raw, p_rel)
#combine pie charts
grid.arrange(nrow=2, pie_all, pie_alls)
#compare relative barplots
grid.arrange(nrow = 2, p_rel, p_rels)



############Visualize PCA for batch/environment effects NOT WORKING###############
# visualize PCA results based on duplicate ID for batch/location/concentration effects
plotIndiv(pca_res, ellipse = TRUE, group = sample_meta_conservative$indiv$concentration, title = "PCA Plot by Concentration", ind.names = clr_data$indiv$dupl_id, legend = TRUE, cex = 0.75, style = "lattice")
clr_data$indiv$indiv$location
plotIndiv(pca.result, group = Data.16S$indiv$concentration, title = "PCA Plot by Concentration", ind.names = Data.16S$indiv$dupl_id, legend = TRUE, cex = 0.75, style = "lattice")
Data.16S$indiv$indiv$location
plot_pca <- function(data.clr, sample.meta, group_var="category", ind_names=TRUE){
  pca.res <- pca(data.clr, logratio="none")
  plotIndiv(pca.res, group=sample.meta[[group_var]], ind.names=ind_names,
            title="PCA Plot (CLR-transformed)", ellipse=TRUE, legend=TRUE)
  return(pca.res)
}
# 3. PCA
pca_res_ext <- plot_pca(clr_data_ext$data.clr, clr_data_ext$sample.meta)


nrow(pca_res$variates$X)     # number of samples in the PCA object
length(sample_meta_conservative$indiv$concentration)
length(sample_meta_conservative$indiv$dupl_id)
sample_meta


rownames(pca_res$variates$X)
meta_aligned <- sample_meta[match(rownames(pca_res$variates$X), sample_meta$SampleID), ]
nrow(meta_aligned)
# should print 22
plotIndiv(
  pca_res,
  ellipse = TRUE,
  group = meta_aligned$concentration,
  ind.names = meta_aligned$submitted.name,   # or dupl_id if you prefer
  title = "PCA Plot by Concentration",
  legend = TRUE,
  cex = 0.75,
  style = "lattice"
)
head(sample_meta)
colnames(sample_meta)
colnames(sample_meta)[1] <- "dupl_id"
sample_meta <- tibble::rownames_to_column(sample_meta, var = "dupl_id")
sample_meta


########Ancom BC
library(ANCOMBC)
library(microbiome)
#variables are ps1_agg and ps1s_agg for conservative and extensive filtering, respectively
# do not use normalized data, enter raw counts
res<- ancombc2(
  data = ps1_agg,
  tax_level = "Genus",  # Specify the taxonomic level
  fix_formula = "category",  # Fixed effects (e.g., treatment groups)
  p_adj_method = "BH",  # P-value adjustment method
  prv_cut = 0.01  # Prevalence cutoff
)
# View all results
head(res$res)

# Filter taxa that are flagged as differentially abundant
sig_taxa <- res$res[res$res$diff_categoryB, ]

# Sort by log fold change
sig_taxa <- sig_taxa[order(sig_taxa$lfc_categoryB, decreasing = TRUE), ]

# Preview top taxa
head(sig_taxa[, c("taxon", "lfc_categoryB", "q_categoryB", "diff_categoryB")])

# Quick summary of significance
table(res$res$diff_categoryB)

# Top 10 taxa by positive effect
top_pos <- head(sig_taxa[order(sig_taxa$beta, decreasing = TRUE), ], 10)

# Top 10 taxa by negative effect
top_neg <- head(sig_taxa[order(sig_taxa$beta, decreasing = FALSE), ], 10)

#Volcano plot
library(ggplot2)
ggplot(res$res, aes(x = lfc_categoryB, y = -log10(q_categoryB))) +
  geom_point(aes(color = diff_categoryB), alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size (log fold change, category B vs A)",
       y = "-log10(adjusted p-value)",
       color = "Significant") +
  theme_minimal() +
  ggtitle("ANCOM-BC2 Volcano Plot")

# Heatmap of significant taxa
library(pheatmap)
# Extract abundance of significant taxa
sig_abund <- t(otu_table(deconps)[rownames(sig_taxa), ])
# Scale abundance for visualization
sig_abund_scaled <- t(scale(t(sig_abund)))

pheatmap(sig_abund_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(sample_data(deconps)),
         main = "Significant Taxa Abundance")


#quick summary statistics
# Number of significant taxa
nrow(sig_taxa)

# Mean effect size of significant taxa
mean(sig_taxa$beta)

# Most positive and negative effect sizes
sig_taxa[which.max(sig_taxa$beta), ]
sig_taxa[which.min(sig_taxa$beta), ]

#write excel with significant taxa
write.csv(sig_taxa, "ANCOMBC2_significant_taxa.csv", row.names = TRUE)


# If NO taxa are signficant (which is true in our case)
# Sort by log fold change regardless of significance
top_taxa <- res$res[order(res$res$lfc_categoryB, decreasing = TRUE), ]
head(top_taxa[, c("taxon", "lfc_categoryB", "q_categoryB", "diff_categoryB")])
# Sort by raw p-value (smallest first)
top_raw_p <- res$res[order(res$res$p_categoryB, decreasing = FALSE), ]
# View top 10
head(top_raw_p[, c("taxon", "lfc_categoryB", "p_categoryB", "q_categoryB", "diff_categoryB")], 10)


####Krona Plot testing#######
# Install psadd from GitHub #could not get krona to work in windows environment
# workaround is to make a krona compatible file to drop into github

library(dplyr)
library(tidyr)
library(readr)
#ensure dupl_id exists
sample_data(ps_filtered)$dupl_id <- rownames(sample_data(ps_filtered))
# Read OTU, taxonomy, and metadata
otu <- as.data.frame(otu_table(ps_filtered))
tax <- as.data.frame(tax_table(ps_filtered))
meta <- as.data.frame(sample_data(ps_filtered))
# If OTUs are in columns, transpose
if(!taxa_are_rows(ps_filtered)){
  otu <- t(otu)
  otu <- as.data.frame(otu)  # <- this ensures it’s a data frame again
}
# Add OTU IDs as a column
otu <- otu %>% 
  rownames_to_column(var = "OTU_ID")
#old: tax <- tax %>%
  #rownames_to_column(var = "OTU_ID") %>%
  #unite("Taxonomy", Kingdom:Species_exact, sep = ";", na.rm = TRUE)
tax <- tax %>%
  rownames_to_column(var = "OTU_ID") %>%
  mutate(across(Kingdom:Species_exact, ~ifelse(is.na(.), "unclassified", .))) %>%
  unite("Taxonomy", Kingdom:Species_exact, sep = ";")
# Combine OTU and taxonomy
otu_tax <- otu %>%
  left_join(tax, by = "OTU_ID")
# Identify which columns in OTU are samples
sample_cols <- intersect(colnames(otu_tax), rownames(meta))
# Pivot longer *only over sample columns*
otu_long <- otu_tax %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "dupl_id",   #adjust to your sample name variable
    values_to = "Count"
  )

# Merge with metadata
otu_meta <- otu_long %>%
  left_join(meta, by = "dupl_id") #adjust to your sample name variable

# Sum counts by group
grouped <- otu_meta %>%
  group_by(category, Taxonomy) %>%
  summarise(Total = sum(Count), .groups = "drop")

# Pivot so each group is a column
krona_input <- grouped %>%
  pivot_wider(names_from = category, values_from = Total, values_fill = 0)

# Write to file for Krona
write_tsv(krona_input, "krona_input_by_group.txt")
write_tsv(dplyr::select(krona_input, A, Taxonomy), "krona_input_by_group.txt")
#####Optional, split file by category to view independently######
# Split by category and write separate files for each
unique_categories <- unique(grouped$category)

for (cat in unique_categories) {
  # Filter for one group
  cat_data <- otu_meta %>%
    filter(category == cat) %>%
    group_by(Taxonomy) %>%
    summarise(Total = sum(Count), .groups = "drop")
  
  # Write to text file
  out_file <- paste0("krona_input_", cat, ".txt")
  write_tsv(cat_data, out_file)
  
  # Optionally generate Krona plot
  cmd <- sprintf('perl "C:/KronaTools/bin/ktImportText.pl" -o "krona_plot_%s.html" "%s"', cat, out_file)
  system(cmd)
}




#testing adding an "other"category
# Read and clean the input
krona_input <- read_tsv("krona_input_A.txt") %>%
  filter(Total > 0) %>%
  mutate(
    Taxonomy = gsub(";+",";", Taxonomy),   # remove repeated semicolons
    Taxonomy = gsub(";+$","", Taxonomy)    # remove trailing semicolons
  )

# Calculate relative abundance
total_counts <- sum(krona_input$Total)
krona_input <- krona_input %>%
  mutate(RelAbund = Total / total_counts * 100)

# Separate "Other" taxa (<1% abundance)
krona_main <- krona_input %>%
  filter(RelAbund >= 1) %>%
  dplyr::select(Total, Taxonomy)

krona_other <- krona_input %>%
  filter(RelAbund < 1) %>%
  summarise(Total = sum(Total)) %>%
  mutate(Taxonomy = "Other") %>%
  dplyr::select(Total, Taxonomy)

# Combine
krona_cleaned <- bind_rows(krona_main, krona_other)

# Write Krona-ready file
write_tsv(dplyr::select(krona_cleaned, Total, Taxonomy), "krona_input_cleaned.txt", col_names = FALSE)

# Generate Krona plot
krona_perl_path <- "C:/KronaTools/bin/ktImportText.pl"
output_html <- "krona_test.html"
cmd <- sprintf('perl "%s" -o "%s" "krona_input_cleaned.txt"', krona_perl_path, output_html)
system(cmd)
cat("Krona plot generated:", output_html, "\n")

#troubleshoot:
read_tsv("krona_input_A.txt", n_max = 10)


# Read and clean the input
krona_input <- read_tsv("krona_input_A.txt") %>%
  filter(Total > 0) %>%
  mutate(
    Taxonomy = gsub(";+",";", Taxonomy),   # remove repeated semicolons
    Taxonomy = gsub(";+$","", Taxonomy)    # remove trailing semicolons
  )
write_tsv(dplyr::select(krona_input,Total, Taxonomy), "krona_input_cleaned.txt", col_names = FALSE)
krona_perl_path <- "C:/KronaTools/bin/ktImportText.pl"
output_html <- "krona_test.html"

cmd <- sprintf('perl "%s" -o "%s" "krona_input_cleaned.txt"', krona_perl_path, output_html)
system(cmd)
cat("Krona plot generated:", output_html, "\n")



#-------------------------------------------
# Intra-Sample Variability Analysis
#-------------------------------------------

#phyloseq object names
  #ps1 = conservative filtering
  #ps1s = extensive filtering
#submitted.name = replicate ID
#dupl_id = sample ID
ps_rel <- transform_sample_counts(ps1s, function(x) x / sum(x))
#aggregate to family level
ps_family <- tax_glom(ps_rel, taxrank = "Family")
#Extract mean abundance across samples
otu_bar <- as(otu_table(ps_family), "matrix") 
if(taxa_are_rows(ps_family)) otu_bar <- t(otu_bar)
#Compute mean abundance and filter by 1%
mean_abund <-colMeans(otu_bar)
keep_taxa <- names(mean_abund[mean_abund > 0.01])
ps1_filt <- prune_taxa(keep_taxa, ps_family)
# Bar plot of relative abundance to see intra-sample variability
plot_bar(ps1_filt, x = "submitted.name", fill = "Family") +
  facet_wrap(~ dupl_id, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Ordination to visualize replicate clustering
# Ordination using Bray-Curtis
ord <- ordinate(ps1_filt, method = "NMDS", distance = "bray")
# Plot ordination colored by sample ID and shaped by replicate ID
plot_ordination(ps1_filt, ord, color = "submitted.name") +
  facet_wrap(~ dupl_id) +
  ggtitle("NMDS Ordination Faceted by Sample") +
  theme_minimal()
plot_ordination(ps1_filt, ord, color = "submitted.name") +
  facet_wrap(~ category) +
  ggtitle("NMDS Ordination Faceted by Sample") +
  theme_minimal()


#Quantify Intra-Sample Dissimilarity
# Extract OTU table
otu <- as(otu_table(ps1_filt), "matrix")
if(taxa_are_rows(ps1_filt)) otu <- t(otu)
# Compute Bray-Curtis dissimilarity
bray <- vegdist(otu, method = "bray")
# Convert to matrix
bray_mat <- as.matrix(bray)
# Convert to long format
library(reshape2)
bray_df_long <- melt(bray_mat, varnames = c("submitted.name.x", "submitted.name.y"), value.name = "bray_dist")
#Prepare metadata
meta <- as(sample_data(ps1_filt), "data.frame")
meta$submitted.name <- rownames(meta)
# Create a second copy for merging the second replicate
meta2 <- meta
meta2$submitted.name.y <- meta2$submitted.name
#Merge metadata
bray_df_long <- merge(bray_df_long, meta, by.x = "submitted.name.x", by.y = "submitted.name")
bray_df_long <- merge(bray_df_long, meta2, by.x = "submitted.name.y", by.y = "submitted.name.y", suffixes = c(".x", ".y"))
# Filter for within-sample comparisons
bray_within <- subset(bray_df_long, dupl_id.x == dupl_id.y & submitted.name.x != submitted.name.y)
# Summarize mean dissimilarity per sample
bray_summary <- bray_within %>%
  group_by(dupl_id.x) %>%
  summarise(mean_dissimilarity = mean(bray_dist))
# Plot
ggplot(bray_summary, aes(x = dupl_id.x, y = mean_dissimilarity)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = round(mean_dissimilarity, 3)), 
           vjust = -0.5, size = 3) +  # Adjust position and size
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Mean Bray-Curtis Dissimilarity") +
  xlab("Sample ID") +
  ggtitle("Intra-Sample Bray-Curtis Dissimilarity")
#bray curtis heatmap
library(pheatmap)
bray <- vegdist(otu, method = "bray")
bray_mat <- as.matrix(bray)
annotation_df <- data.frame(dupl_id = sample_data(ps1_filt)$dupl_id)
rownames(annotation_df) <- rownames(sample_data(ps1_filt))
pheatmap(bray_mat, annotation_col = annotation_df)

# Alpha-diversity per replicate
alpha_div <- estimate_richness(ps1, measures = c("Shannon", "Simpson"))
alpha_div$sample_id <- rownames(alpha_div)
# Prepare metadata
meta <- as(sample_data(ps1), "data.frame")
meta$sample_id <- rownames(meta)  # Also P30859_101 etc.
# Merge with metadata
alpha_div <- merge(alpha_div, meta, by = "sample_id")
alpha_div$dupl_id <- as.factor(alpha_div$dupl_id)
# Plot Shannon diversity grouped by sample
ggplot(alpha_div, aes(x = dupl_id, y = Shannon)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Shannon Diversity") +
  xlab("Sample ID") +
  ggtitle("Alpha Diversity per Replicate")


# Bar plot by location on aggregated samples
sample_data(ps1_agg)$dupl_id <- rownames(sample_data(ps1_agg))
ps_rel_agg <- transform_sample_counts(ps1_agg, function(x) x / sum(x))
ps_family_agg <- tax_glom(ps_rel_agg, taxrank = "Family")
otu_bar_agg <- as(otu_table(ps_family_agg), "matrix")
if(taxa_are_rows(ps_family_agg)) otu_bar_agg <- t(otu_bar_agg)
mean_abund_agg <-colMeans(otu_bar_agg)
keep_taxa <- names(mean_abund_agg[mean_abund_agg > 0.01])
ps1_agg_filt <- prune_taxa(keep_taxa, ps_family_agg)
plot_bar(ps1_agg_filt, x = "dupl_id", fill = "Family") +
  facet_wrap(~ location, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#-------------Intra-sample variability as a function--------

run_intra_sample_analysis <- function(ps_object, object_name = "ps") {
  # Relative abundance
  ps_rel <- transform_sample_counts(ps_object, function(x) x / sum(x))
  
  # Aggregate to Family level
  ps_family <- tax_glom(ps_rel, taxrank = "Family")
  
  # Filter by mean abundance > 1%
  otu_bar <- as(otu_table(ps_family), "matrix")
  if (taxa_are_rows(ps_family)) otu_bar <- t(otu_bar)
  mean_abund <- colMeans(otu_bar)
  keep_taxa <- names(mean_abund[mean_abund > 0.01])
  ps_filt <- prune_taxa(keep_taxa, ps_family)
  
  # Bar plot
  print(
    plot_bar(ps_filt, x = "submitted.name", fill = "Family") +
      facet_wrap(~ dupl_id, scales = "free_x") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(paste("Relative Abundance -", object_name))
  )
  
  # Ordination
  ord <- ordinate(ps_filt, method = "NMDS", distance = "bray")
  print(
    plot_ordination(ps_filt, ord, color = "submitted.name") +
      facet_wrap(~ dupl_id) +
      ggtitle(paste("NMDS Ordination by dupl_id -", object_name)) +
      theme_minimal()
  )
  print(
    plot_ordination(ps_filt, ord, color = "submitted.name") +
      facet_wrap(~ category) +
      ggtitle(paste("NMDS Ordination by category -", object_name)) +
      theme_minimal()
  )
  
  # Bray-Curtis dissimilarity
  otu <- as(otu_table(ps_filt), "matrix")
  if (taxa_are_rows(ps_filt)) otu <- t(otu)
  bray <- vegdist(otu, method = "bray")
  bray_mat <- as.matrix(bray)
  
  # Metadata
  meta <- as(sample_data(ps_filt), "data.frame")
  meta$submitted.name <- rownames(meta)
  meta2 <- meta
  meta2$submitted.name.y <- meta2$submitted.name
  
  # Merge and filter
  bray_df_long <- reshape2::melt(bray_mat, varnames = c("submitted.name.x", "submitted.name.y"), value.name = "bray_dist")
  bray_df_long <- merge(bray_df_long, meta, by.x = "submitted.name.x", by.y = "submitted.name")
  bray_df_long <- merge(bray_df_long, meta2, by.x = "submitted.name.y", by.y = "submitted.name.y", suffixes = c(".x", ".y"))
  bray_within <- subset(bray_df_long, dupl_id.x == dupl_id.y & submitted.name.x != submitted.name.y)
  
  # Summary
  bray_summary <- bray_within %>%
    group_by(dupl_id.x) %>%
    summarise(mean_dissimilarity = mean(bray_dist))
  
  # Plot
  print(
    ggplot(bray_summary, aes(x = dupl_id.x, y = mean_dissimilarity)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_text(aes(label = round(mean_dissimilarity, 3)), vjust = -0.5, size = 3) +
      theme(axis.text.x = element_text(angle = 90)) +
      ylab("Mean Bray-Curtis Dissimilarity") +
      xlab("Sample ID") +
      ggtitle(paste("Intra-Sample Dissimilarity -", object_name))
  )
  
  # Heatmap
  annotation_df <- data.frame(dupl_id = sample_data(ps_filt)$dupl_id)
  rownames(annotation_df) <- rownames(sample_data(ps_filt))
  pheatmap(bray_mat, annotation_col = annotation_df)
  
  # Alpha diversity
  alpha_div <- estimate_richness(ps_object, measures = c("Shannon", "Simpson"))
  alpha_div$sample_id <- rownames(alpha_div)
  meta_alpha <- as(sample_data(ps_object), "data.frame")
  meta_alpha$sample_id <- rownames(meta_alpha)
  alpha_div <- merge(alpha_div, meta_alpha, by = "sample_id")
  alpha_div$dupl_id <- as.factor(alpha_div$dupl_id)
  
  # Plot alpha diversity
  print(
    ggplot(alpha_div, aes(x = dupl_id, y = Shannon)) +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 90)) +
      ylab("Shannon Diversity") +
      xlab("Sample ID") +
      ggtitle(paste("Alpha Diversity per Replicate -", object_name))
  )
}

run_intra_sample_analysis(ps1, "Conservative Filtering")

#troubshooting ps1s, which isn't running as expected
# Clean Genus first
tax_tab <- tax_table(ps1)
genus <- as.character(tax_tab[, "Genus"])
genus[is.na(genus) | genus == ""] <- "Unclassified_Genus"
tax_table(ps1) <- tax_table(matrix(genus, ncol=1, dimnames=list(rownames(tax_tab), "Genus")))

# Filter all at once
ps1s <- subset_taxa(ps1, !(Genus %in% c(filter_genera, filter_contaminant_taxa)))

run_intra_sample_analysis(ps1s, "Extensive Filtering")




#-------------------------
#Other intra-sample options
#--------------------------
library(dplyr)

# Calculate summary stats per sample group
alpha_summary <- alpha_div %>%
  group_by(dupl_id) %>%
  summarise(
    mean_shannon = mean(Shannon),
    sd_shannon = sd(Shannon),
    cv_shannon = sd_shannon / mean_shannon,
    mean_simpson = mean(Simpson),
    sd_simpson = sd(Simpson),
    cv_simpson = sd_simpson / mean_simpson
  )
#Detect outliers
# z-score
alpha_div <- alpha_div %>%
  group_by(dupl_id) %>%
  mutate(
    z_shannon = ifelse(sd(Shannon) > 0, (Shannon - mean(Shannon)) / sd(Shannon), 0),
    z_simpson = ifelse(sd(Simpson) > 0, (Simpson - mean(Simpson)) / sd(Simpson), 0)
  )

z_outliers <- alpha_div %>%
  filter(abs(z_shannon) > 2)
# IQR method---this is better for our purposes
iqr_outliers <- alpha_div %>%
  group_by(dupl_id) %>%
  filter(
    Shannon < quantile(Shannon, 0.25) - 1.5 * IQR(Shannon) |
      Shannon > quantile(Shannon, 0.75) + 1.5 * IQR(Shannon)
  )
#z-score outlier plot
ggplot(alpha_div, aes(x = dupl_id, y = Shannon)) +
  geom_boxplot(outlier.shape = NA) +  # Hide default outliers
  geom_jitter(width = 0.2, alpha = 0.6) +  # Show all points
  geom_point(data = z_outliers, aes(x = dupl_id, y = Shannon), color = "red", size = 3) +
  geom_text(data = z_outliers, aes(label = sample_id), vjust = -0.5, size = 2.5) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Alpha Diversity with Z-score Outliers Highlighted")
#iqr outlier plot
ggplot(alpha_div, aes(x = dupl_id, y = Shannon)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_point(data = iqr_outliers, aes(x = dupl_id, y = Shannon), color = "blue", size = 3) +
  geom_text(data = iqr_outliers, aes(label = sample_id), vjust = -0.5, size = 2.5) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Alpha Diversity with IQR Outliers Highlighted")

#Boxplot with outlier labels
ggplot(alpha_div, aes(x = dupl_id, y = Shannon)) +
  geom_boxplot(outlier.colour = "red", outlier.shape = 8) +
  geom_text(data = iqr_outliers, aes(label = sample_id), vjust = -0.5, size = 2.5) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Shannon Diversity") +
  ggtitle("Alpha Diversity with Outliers")

#beta diveristy outliers


bray_summary <- bray_within %>%
  group_by(dupl_id.x) %>%
  summarise(mean_dissimilarity = mean(bray_dist),
            sd_dissimilarity = sd(bray_dist))

bray_within <- bray_within %>%
  group_by(dupl_id.x) %>%
  mutate(rel_dev = abs(bray_dist - mean(bray_dist)) / mean(bray_dist))
bray_outliers <- bray_within %>% filter(rel_dev > 0.5)

ggplot(bray_within, aes(x = dupl_id.x, y = bray_dist)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_point(data = bray_outliers, aes(x = dupl_id.x, y = bray_dist), color = "red", size = 3) +
  geom_text(data = bray_outliers, aes(label = paste(submitted.name.x, submitted.name.y, sep = "-")), 
            vjust = -0.5, size = 2.5) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Bray-Curtis Dissimilarity") +
  ggtitle("Replicate Pair Dissimilarity with Outliers Highlighted")

ggplot(bray_summary, aes(x = dupl_id.x, y = mean_dissimilarity, fill = mean_dissimilarity > 0.6)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
  geom_text(aes(label = round(mean_dissimilarity, 3)), vjust = -0.5, size = 3) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Mean Bray-Curtis Dissimilarity") +
  xlab("Sample ID") +
  ggtitle("Highlight Samples with Mean Bray-Curtis > 0.6") +
  guides(fill = FALSE)  # Hide legend if not needed

######Things to do still
#  1) Intra-sample variability analysis
#  2) Deciding sPLS-DA prevalence/abundance filtering for pilot study
#  3) Finalize ANCOM-BC2 analysis
#  4) Make Krona Plot (in testing)