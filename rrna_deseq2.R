#!/usr/bin/env Rscript

required_packages <- c("DESeq2", "data.table", "ggplot2", "pheatmap", "RColorBrewer", "dplyr", "tidyr", "ggrepel", "VennDiagram", "grid")
new_packages <- required_packages[!(required_packages %in% rownames(installed.packages()))]
if(length(new_packages) > 0) install.packages(new_packages, repos = "https://cloud.r-project.org")
invisible(lapply(required_packages, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))))

args <- commandArgs(trailingOnly = TRUE)
otu_file <- args[1]
metadata_file <- args[2]
output_dir <- args[3]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

otu_df <- fread(otu_file)
otu_df <- as.data.frame(otu_df)
rownames(otu_df) <- otu_df[[1]]
otu_df[[1]] <- NULL

metadata <- fread(metadata_file)
metadata <- as.data.frame(metadata)
names(metadata) <- tolower(names(metadata))
sample_col <- grep("sample", names(metadata), value = TRUE)[1]
cond_col   <- grep("condition", names(metadata), value = TRUE)[1]
colnames(metadata)[colnames(metadata) == sample_col] <- "SampleID"
colnames(metadata)[colnames(metadata) == cond_col] <- "Condition"

total_counts <- rowSums(otu_df)
prevalence <- rowSums(otu_df > 0) / ncol(otu_df)
filtered_otu <- otu_df[total_counts >= 20 & prevalence >= 0.10, , drop = FALSE]

samples_in_both <- intersect(colnames(filtered_otu), metadata$SampleID)
filtered_otu <- filtered_otu[, samples_in_both, drop = FALSE]
metadata <- metadata[metadata$SampleID %in% samples_in_both, , drop = FALSE]
metadata <- metadata[match(samples_in_both, metadata$SampleID), , drop = FALSE]

metadata$Condition <- factor(metadata$Condition, levels = c("NoHediste_NoNoise", "Hediste_NoNoise", "NoHediste_Noise", "Hediste_Noise"))

dds <- DESeqDataSetFromMatrix(countData = round(filtered_otu),
                              colData = metadata,
                              design = ~ Condition)

dds <- DESeq(dds)
rld <- rlog(dds, blind = TRUE)

contrast_list <- list(
  Hediste_vs_NoHediste_NoNoise = c("Condition", "Hediste_NoNoise", "NoHediste_NoNoise"),
  Noise_vs_NoNoise_NoHediste = c("Condition", "NoHediste_Noise", "NoHediste_NoNoise"),
  Noise_vs_NoNoise_Hediste = c("Condition", "Hediste_Noise", "Hediste_NoNoise")
)

sig_genus_sets <- list()

for (contrast_name in names(contrast_list)) {
  contrast <- contrast_list[[contrast_name]]
  res <- results(dds, contrast = contrast)
  resOrdered <- res[order(res$padj), ]
  res_df <- as.data.frame(resOrdered)
  res_df$Genus <- rownames(res_df)
  fwrite(res_df, file = file.path(output_dir, paste0("deseq2_results_", contrast_name, ".tsv")), sep = "\t")

  png(filename = file.path(output_dir, paste0("MA_plot_", contrast_name, ".png")), width = 800, height = 600)
  plotMA(res, main = paste("DESeq2 MA-plot:", contrast_name), ylim = c(-5, 5))
  dev.off()

  sig_res <- res_df %>% filter(!is.na(padj) & padj <= 0.05)
  sig_genus_sets[[contrast_name]] <- sig_res$Genus

  upregulated <- sig_res %>% filter(log2FoldChange > 1)
  downregulated <- sig_res %>% filter(log2FoldChange < -1)

  fwrite(upregulated, file = file.path(output_dir, paste0("upregulated_genera_", contrast_name, ".tsv")), sep = "\t")
  fwrite(downregulated, file = file.path(output_dir, paste0("downregulated_genera_", contrast_name, ".tsv")), sep = "\t")

  annotation_col <- as.data.frame(colData(rld)[ , "Condition", drop = FALSE])

  if (nrow(upregulated) > 1) {
    top_up <- upregulated %>% arrange(desc(log2FoldChange)) %>% head(25) %>% pull(Genus)
    mat_up <- assay(rld)[top_up, , drop = FALSE]
    if (nrow(mat_up) >= 2) {
      mat_up <- mat_up - rowMeans(mat_up)
      png(filename = file.path(output_dir, paste0("heatmap_upregulated_", contrast_name, ".png")), width = 1000, height = 800)
      pheatmap(mat_up,
               annotation_col = annotation_col,
               fontsize_row = 10,
               fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               show_rownames = TRUE,
               show_colnames = TRUE,
               main = paste("Top Upregulated Genera:", contrast_name))
      dev.off()
    }
  }

  if (nrow(downregulated) > 1) {
    top_down <- downregulated %>% arrange(log2FoldChange) %>% head(25) %>% pull(Genus)
    mat_down <- assay(rld)[top_down, , drop = FALSE]
    if (nrow(mat_down) >= 2) {
      mat_down <- mat_down - rowMeans(mat_down)
      png(filename = file.path(output_dir, paste0("heatmap_downregulated_", contrast_name, ".png")), width = 1000, height = 800)
      pheatmap(mat_down,
               annotation_col = annotation_col,
               fontsize_row = 10,
               fontsize_col = 10,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               clustering_method = "complete",
               show_rownames = TRUE,
               show_colnames = TRUE,
               main = paste("Top Downregulated Genera:", contrast_name))
      dev.off()
    }
  }

  png(filename = file.path(output_dir, paste0("volcano_plot_", contrast_name, ".png")), width = 800, height = 600)
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = padj <= 0.05)) +
    geom_text_repel(data = subset(res_df, padj <= 0.05 & abs(log2FoldChange) > 2), aes(label = Genus), size = 3) +
    scale_color_manual(values = c("gray", "red")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", contrast_name), x = "log2 Fold Change", y = "-log10 Adjusted P-value")
  dev.off()
}

venn_out <- file.path(output_dir, "venn_significant_genera.png")
png(venn_out, width = 800, height = 600)
venn.plot <- venn.diagram(
  sig_genus_sets,
  filename = NULL,
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("red", "blue", "green"),
  main = "Significant Genera Overlap"
)
grid.draw(venn.plot)
dev.off()

rel_abund <- sweep(filtered_otu, 2, colSums(filtered_otu), FUN = "/")
top_genera_bar <- names(sort(rowSums(rel_abund), decreasing = TRUE))[1:20]
rel_abund_top <- rel_abund[top_genera_bar, , drop = FALSE]
rel_abund_top_df <- as.data.frame(t(rel_abund_top))
rel_abund_top_df$SampleID <- rownames(rel_abund_top_df)
rel_abund_long <- pivot_longer(rel_abund_top_df, cols = -SampleID, names_to = "Genus", values_to = "RelativeAbundance")
rel_abund_long <- merge(rel_abund_long, metadata[, c("SampleID", "Condition")], by = "SampleID")

png(filename = file.path(output_dir, "stacked_barplot_top_genera.png"), width = 1000, height = 600)
ggplot(rel_abund_long, aes(x = SampleID, y = RelativeAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Condition, scales = "free_x") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Relative Abundance of Top Genera", y = "Relative Abundance", x = "Sample") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(20))
dev.off()

rel_abund_avg <- rel_abund_long %>%
  group_by(Condition, Genus) %>%
  summarise(MeanRelAbund = mean(RelativeAbundance), .groups = "drop")

png(filename = file.path(output_dir, "stacked_barplot_top_genera_by_condition.png"), width = 800, height = 600)
ggplot(rel_abund_avg, aes(x = Condition, y = MeanRelAbund, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 12) +
  labs(title = "Mean Relative Abundance of Top Genera by Condition", y = "Relative Abundance", x = "Condition") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(20))
dev.off()
