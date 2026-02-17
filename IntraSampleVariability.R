library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(reshape2)
run_intra_sample_analysis <- function(ps_object, object_name = "ps", taxrank = "Family",
                                      abundance_cutoff = 0.01,
                                      absolute_bc_cutoff = 0.6,
                                      ord_method = c("PCoA", "NMDS")) {
  
  ord_method <- match.arg(ord_method)
  
  ## -------------------------
  ## 0. Defensive cleanup
  ## -------------------------
  ps_object <- prune_taxa(taxa_sums(ps_object) > 0, ps_object)
  ps_object <- prune_samples(sample_sums(ps_object) > 0, ps_object)
  
  ## -------------------------
  ## 0.5 Ensure submitted.name exists
  ## -------------------------
  sd <- as.data.frame(sample_data(ps_object), check.names = FALSE)
  
  if (!"submitted.name" %in% colnames(sd)) {
    sd$submitted.name <- rownames(sd)
    sample_data(ps_object) <- sample_data(sd)
  }
  
  ## -------------------------
  ## 1. Aggregate to Family
  ## -------------------------
  ps_fam <- tax_glom(ps_object, taxrank = taxrank)
  
  ## Relative abundance
  ps_rel <- transform_sample_counts(ps_fam, function(x) x / sum(x))
  
  ## -------------------------
  ## 2. Filter low-abundance taxa
  ## -------------------------
  otu <- as(otu_table(ps_rel), "matrix")
  if (taxa_are_rows(ps_rel)) otu <- t(otu)
  
  mean_abund <- colMeans(otu)
  keep_taxa <- names(mean_abund[mean_abund > abundance_cutoff])
  
  # Conditional pruning
  if(length(keep_taxa) == ntaxa(ps_rel)) {
    # If all taxa are kept, skip prune_taxa
    ps_filt <- ps_rel
  } else {
    # Otherwise, prune normally
    ps_filt <- prune_taxa(keep_taxa, ps_rel)
  }
  
  ps_filt <- prune_samples(sample_sums(ps_filt) > 0, ps_filt)
  
  # Stop Message
  if (ntaxa(ps_filt) == 0 | nsamples(ps_filt) == 0) {
    stop("No taxa or samples remaining after filtering. Consider lowering abundance_cutoff or using less aggressive filtering.")
  }
  # 
  
  
  # -------------------------
  # 2.5 Bar plot of relative abundance
  # -------------------------
  
  print(
    plot_bar(ps_filt, x = "submitted.name", fill = taxrank) +
      facet_wrap(~ dupl_id, scales = "free_x") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none"
      ) +
      ggtitle(paste("Relative Abundance -", object_name))
  )
  
  
  # -------------------------
  # 3. Ordination (PCoA or NMDS)
  # -------------------------
  otu <- as(otu_table(ps_filt), "matrix")
  if(taxa_are_rows(ps_filt)) otu <- t(otu)
  
  if(ord_method == "PCoA") {
    ord <- cmdscale(vegdist(otu, method = "bray"), k = 2, eig = TRUE)
    ord_df <- data.frame(ord1 = ord$points[,1], ord2 = ord$points[,2], dupl_id = sample_data(ps_filt)$dupl_id)
  } else {
    ord <- metaMDS(otu, distance = "bray", k = 2, trymax = 50)
    ord_df <- data.frame(ord1 = ord$points[,1], ord2 = ord$points[,2], dupl_id = sample_data(ps_filt)$dupl_id)
  }
  
  print(
    ggplot(ord_df, aes(x = ord1, y = ord2, color = dupl_id)) +
      geom_point(size = 3) +
      geom_text(aes(label = dupl_id), vjust = -0.5, size = 5) +  # Add labels
      ggtitle(paste("Ordination (", ord_method, ") - replicates colored by sample", sep = "")) +
      theme_minimal()
  )
  
  # -------------------------
  # 4. Intra-sample Bray-Curtis dissimilarity
  # -------------------------
  bray <- vegdist(otu, method = "bray")
  bray_mat <- as.matrix(bray)
  
  # Metadata
  meta <- as(sample_data(ps_filt), "data.frame")
  meta$submitted.name <- rownames(meta)
  meta2 <- meta
  meta2$submitted.name.y <- meta2$submitted.name
  
  # Long format + within-sample comparisons
  bray_df_long <- reshape2::melt(bray_mat, varnames = c("submitted.name.x", "submitted.name.y"), value.name = "bray_dist")
  bray_df_long <- merge(bray_df_long, meta, by.x = "submitted.name.x", by.y = "submitted.name")
  bray_df_long <- merge(bray_df_long, meta2, by.x = "submitted.name.y", by.y = "submitted.name.y", suffixes = c(".x", ".y"))
  bray_within <- subset(bray_df_long, dupl_id.x == dupl_id.y & submitted.name.x != submitted.name.y)
  
  bray_summary <- bray_within %>%
    group_by(dupl_id.x) %>%
    summarise(mean_dissimilarity = mean(bray_dist),
              sd_dissimilarity = sd(bray_dist))
  
  # -------------------------
  # 5. Highlight samples exceeding absolute BC cutoff
  # -------------------------
  print(
    ggplot(bray_summary, aes(x = dupl_id.x, y = mean_dissimilarity, fill = mean_dissimilarity > absolute_bc_cutoff)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "red")) +
      geom_text(aes(label = round(mean_dissimilarity, 3)), vjust = -0.5, size = 3) +
      theme(axis.text.x = element_text(angle = 90)) +
      ylab("Mean Bray-Curtis Dissimilarity") +
      xlab("Sample ID") +
      ggtitle(paste("Samples Exceeding Absolute BC Cutoff (>", absolute_bc_cutoff, ")", sep = " ")) +
      guides(fill = FALSE)
  )
  
  # -------------------------
  # 6. Bray-Curtis heatmap
  # -------------------------
  annotation_df <- data.frame(dupl_id = sample_data(ps_filt)$dupl_id)
  rownames(annotation_df) <- rownames(sample_data(ps_filt))
  pheatmap(bray_mat, annotation_col = annotation_df, main = paste("BC Heatmap -", object_name))
  
  # -------------------------
  # 7. Alpha diversity
  # -------------------------
  alpha_div <- estimate_richness(ps_object, measures = c("Shannon", "Simpson"))
  alpha_div$sample_id <- rownames(alpha_div)
  meta_alpha <- as(sample_data(ps_object), "data.frame")
  meta_alpha$sample_id <- rownames(meta_alpha)
  alpha_div <- merge(alpha_div, meta_alpha, by = "sample_id")
  
  print(
    ggplot(alpha_div, aes(x = factor(dupl_id), y = Shannon)) +
      geom_boxplot() +
      theme(axis.text.x = element_text(angle = 90)) +
      ylab("Shannon Diversity") +
      xlab("Sample ID") +
      ggtitle(paste("Alpha Diversity per Replicate -", object_name))
  )
  
  # -------------------------
  # Return results
  # -------------------------
  return(list(
    filtered_phyloseq = ps_filt,
    bray_summary = bray_summary,
    bray_within = bray_within,
    alpha_div = alpha_div,
    ordination_df = ord_df
  ))
}

#usage: 
# Conservative filtering 
res_ps1 <- run_intra_sample_analysis(
  ps_object = ps1,
  object_name = "Conservative Filtering",
  taxrank = "Family",
  abundance_cutoff = 0.01, #change to zero for prettier plot but do not use for BC distances
  absolute_bc_cutoff = 0.5,
  ord_method = "NMDS"  #can change to PCoA
)
names(res_ps1)


# Extensive filtering 
res_ps1s <- run_intra_sample_analysis(
  ps_object = ps1s,
  object_name = "Extensive Filtering",
  taxrank = "Family",
  abundance_cutoff = 0.01,
  absolute_bc_cutoff = 0.5,
  ord_method = "NMDS"  #remove this line if you want PCoA; would recommend nmds generally
)



#investigate alpha diversity:
sample_sums(ps_filt)[ps_filt$dupl_id == "X"]  # for that dupl_id



get_alpha_diff_outliers <- function(res_obj, diff_threshold = 0.5) {
  
  res_obj$alpha_div %>%
    group_by(dupl_id) %>%
    filter(n() == 2) %>%                     # if exactly 2 replicates
    arrange(submitted.name) %>%              # enforce stable order
    reframe(
      shannon_diff = abs(Shannon[1] - Shannon[2]),
      submitted_names = paste(submitted.name, collapse = ", "),
      Shannon_values = paste(round(Shannon, 3), collapse = ", ")
    ) %>%
    filter(shannon_diff > diff_threshold)
}

# Conservative
alpha_diff_outliers_ps1 <- get_alpha_diff_outliers_pct(res_ps1, diff_threshold = 0.5)
alpha_diff_outliers_ps1

# Extensive
alpha_diff_outliers_ps1s <- get_alpha_diff_outliers(res_ps1s, diff_threshold = 0.5)
alpha_diff_outliers_ps1s
