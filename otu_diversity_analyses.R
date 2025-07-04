#!/usr/bin/env Rscript

# ===============================
# Alpha & Beta Diversity Analysis
# ===============================

# -------------------------------
# Install and load required packages
# -------------------------------
packages <- c("vegan", "phyloseq", "ggplot2", "reshape2", "RColorBrewer",
              "gridExtra", "FSA", "rstatix", "ggpubr", "dunn.test",
              "purrr", "tidyr", "dplyr", "data.table")

installed <- rownames(installed.packages())
lapply(setdiff(packages, installed), install.packages, repos = "https://cloud.r-project.org")
invisible(lapply(packages, library, character.only = TRUE))

# -------------------------------
# Parse command-line arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
otu_file <- args[1]
metadata_file <- args[2]
output_dir <- args[3]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# Load OTU Table
# -------------------------------
otu_df <- read.table(otu_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)

# Apply genus-level filtering: total count ≥ 20 and prevalence ≥ 10%
otu_counts <- rowSums(otu_df)
otu_prevalence <- rowSums(otu_df > 0)
sample_threshold <- ceiling(0.10 * ncol(otu_df))
otu_df_filtered <- otu_df[otu_counts >= 20 & otu_prevalence >= sample_threshold, ]

# -------------------------------
# Prepare OTU matrix
# -------------------------------
otu_matrix <- as.matrix(t(otu_df_filtered))
otu_matrix <- otu_matrix[rowSums(otu_matrix) > 0, ]

# Calculate total abundance before rarefaction
TotalAbundance <- rowSums(otu_matrix)

# -------------------------------
# Rarefy to minimum depth
# -------------------------------
min_depth <- min(rowSums(otu_matrix))
set.seed(42)
otu_rare <- rrarefy(otu_matrix, sample = min_depth)

# -------------------------------
# Load and clean metadata
# -------------------------------
metadata <- read.table(metadata_file, sep = "\t", header = TRUE, check.names = FALSE)
names(metadata) <- trimws(tolower(names(metadata)))  # clean and lowercase column names

# Detect columns
sample_col <- grep("sample", names(metadata), value = TRUE)[1]
cond_col   <- grep("condition", names(metadata), value = TRUE)[1]

if (is.na(sample_col) || is.na(cond_col)) {
  stop("\u274c ERROR: Could not detect sample ID or condition columns in metadata.")
}

# -------------------------------
# Alpha Diversity Calculation
# -------------------------------
otu_sample_ids <- rownames(otu_rare)


# Alpha diversity metrics
Observed <- specnumber(otu_rare)
Shannon <- diversity(otu_rare, index = "shannon")
Simpson <- diversity(otu_rare, index = "simpson")  # Simpson index
PielouEvenness <- Shannon / log(Observed)
Chao1 <- estimateR(otu_rare)["S.chao1", ]

# Combine into dataframe
alpha_div <- data.frame(
  SampleID = otu_sample_ids,
  TotalAbundance = TotalAbundance,
  Observed = Observed,
  Chao1 = Chao1,
  Shannon = Shannon,
  Simpson = Simpson,
  PielouEvenness = PielouEvenness
)

# Merge with metadata
metadata_subset <- metadata[, c(sample_col, cond_col)]
colnames(metadata_subset) <- c("SampleID", "Condition")
alpha_div <- merge(alpha_div, metadata_subset, by = "SampleID", all.x = TRUE)

# -------------------------------
# Save Alpha Diversity Table
# -------------------------------
write.table(alpha_div, file = file.path(output_dir, "alpha_diversity_with_conditions.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# -------------------------------
# Summary Stats Table (mean, sd, sem)
# -------------------------------
alpha_long <- alpha_div %>%
  tidyr::pivot_longer(cols = c(Observed, Chao1, Shannon, PielouEvenness),
                      names_to = "Metric", values_to = "Value")

summary_stats <- alpha_long %>%
  group_by(Condition, Metric) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    SEM = SD / sqrt(n()),
    .groups = "drop"
  )

fwrite(summary_stats, file = file.path(output_dir, "alpha_diversity_summary_stats.tsv"), sep = "\t", quote = FALSE)

# -------------------------------
# Plot Alpha Diversity (white background)
# -------------------------------
alpha_melt <- reshape2::melt(alpha_div, id.vars = c("SampleID", "Condition"),
                             variable.name = "Metric", value.name = "Value")

p_alpha <- ggplot(alpha_melt, aes(x = Condition, y = Value, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Alpha Diversity", x = "", y = "Diversity Index") +
  scale_fill_brewer(palette = "Set2")

ggsave(filename = file.path(output_dir, "alpha_diversity_by_condition.png"),
       plot = p_alpha, width = 10, height = 6)

# -------------------------------
# Alpha Diversity Stats
# -------------------------------
alpha_nested <- alpha_melt %>%
  group_by(Metric) %>%
  nest()

alpha_stats <- alpha_nested %>%
  mutate(test = map(data, ~ tryCatch(
    kruskal_test(Value ~ Condition, data = .x),
    error = function(e) data.frame(statistic = NA, df = NA, p = NA, method = "kruskal_test", warning = e$message)
  ))) %>%
  select(Metric, test) %>%
  unnest(cols = test)

dunn_results <- alpha_nested %>%
  mutate(result = map(data, ~ tryCatch(
    dunn_test(Value ~ Condition, data = .x, p.adjust.method = "fdr"),
    error = function(e) data.frame(group1 = NA, group2 = NA, p.adj = NA, p.adj.signif = NA, warning = e$message)
  ))) %>%
  select(Metric, result) %>%
  unnest(cols = result)

write.table(alpha_stats, file = file.path(output_dir, "alpha_kruskal_test.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(dunn_results, file = file.path(output_dir, "alpha_dunn_test.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# -------------------------------
# Beta Diversity (Bray-Curtis PCoA)
# -------------------------------
bray_dist <- vegdist(otu_rare, method = "bray")
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
pcoa_df <- data.frame(SampleID = rownames(pcoa$points),
                      Axis1 = pcoa$points[, 1],
                      Axis2 = pcoa$points[, 2])

# Merge metadata to PCoA
pcoa_df <- merge(pcoa_df, metadata_subset, by = "SampleID", all.x = TRUE)

# Save PCoA coordinates
write.table(pcoa_df, file = file.path(output_dir, "beta_diversity_pcoa_coordinates.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Plot PCoA (white background)
p_beta <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Condition)) +
  stat_ellipse(type = "norm", linetype = 2, size = 1, alpha = 0.5) +
  geom_point(size = 3, alpha = 0.9) +
  theme_classic(base_size = 14) +
  labs(title = "Beta Diversity (Bray-Curtis PCoA)",
       x = "PCoA1", y = "PCoA2") +
  scale_color_brewer(palette = "Set1")

ggsave(filename = file.path(output_dir, "beta_diversity_pcoa.png"),
       plot = p_beta, width = 8, height = 6)

# -------------------------------
# PERMANOVA
# -------------------------------
adonis_result <- tryCatch(
  adonis2(bray_dist ~ Condition, data = pcoa_df),
  error = function(e) paste("PERMANOVA failed:", e$message)
)

capture.output(adonis_result, file = file.path(output_dir, "beta_diversity_permanova.txt"))
