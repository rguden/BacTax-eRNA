#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
out_species <- args[2]
out_genus <- args[3]
threads <- as.integer(args[4])

library(data.table)
library(stringr)
library(future.apply)

plan(multisession, workers = threads)

# Get all Kraken filtered Bacteria report files
kraken_files <- list.files(input_dir, pattern = "_report_bacteria_only\\.txt$", full.names = TRUE)

# Parse all files in parallel
otu_list <- future_lapply(kraken_files, function(file) {
  sample <- str_replace(basename(file), "_report_bacteria_only\\.txt", "")
  df <- fread(file, header = FALSE, sep = "\t", quote = "")
  setnames(df, c("percent", "reads_clade", "reads_taxon", "rank", "taxid", "name"))
  df$name <- trimws(df$name)

  df[, .(sample = sample, rank = rank, taxon = name, count = reads_taxon)]
})

# Combine all samples
otu_all <- rbindlist(otu_list)

# Function to make OTU table for a specific rank
make_otu_table <- function(rank_code, taxon_label, out_path) {
  otu_rank <- otu_all[rank == rank_code]
  otu_wide <- dcast(otu_rank, taxon ~ sample, value.var = "count", fun.aggregate = sum, fill = 0)
  setnames(otu_wide, "taxon", taxon_label)
  fwrite(otu_wide, file = out_path, sep = "\t", quote = FALSE)
}

# Generate Genus and Species OTU tables
make_otu_table("G", "Genus", out_genus)
make_otu_table("S", "Species", out_species)
