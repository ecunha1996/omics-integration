run_gage_analysis <- function(getmm, gmt_path, output_path, control=NULL, sample=NULL, pval_threshold=0.05, bidirectional=T) {
  # Load required libraries
  if (!requireNamespace("gage", quietly = TRUE)) install.packages("gage")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("gage", quietly = TRUE)) BiocManager::install("gage")
  
  library(gage)
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  
  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
  
  cn <- colnames(getmm)
  
  # Read custom pathways
  custom_pathways <- read.gmt(gmt_path)
  
  # Convert pathways to named list
  gsets_list <- custom_pathways %>%
    group_by(term) %>%
    summarise(genes = list(gene)) %>%
    deframe()
  
  # Identify control and treatment samples
  control_3 <- grep(control, cn, ignore.case = TRUE)
  nitrogen_3 <- grep(sample, cn, ignore.case = TRUE)
  
  
  gage_results <- gage(getmm, gsets = gsets_list, ref = control_3, samp = nitrogen_3,
                       compare = "unpaired", same.dir = bidirectional, saaTest = gs.KSTest, set.size=c(15, 50000))
  
  # Extract significant gene sets
  
  sig_results <- sigGeneSet(gage_results, outname = file.path(output_path, paste0("sig_results_", control, "_vs_", sample)))
  
  if (bidirectional){
  
  greater_file <- file.path(output_path, paste0("greater_", control, "_vs_", sample, ".tsv"))
  less_file <- file.path(output_path, paste0("less_", control, "_vs_", sample, ".tsv"))
  
  write.table(gage_results$greater, file = greater_file, sep = "\t", quote = FALSE, row.names = TRUE)
  write.table(gage_results$less, file = less_file, sep = "\t", quote = FALSE, row.names = TRUE)
  
  greater_significant <- gage_results$greater %>% as.data.frame() %>% filter(q.val < pval_threshold)
  less_significant <- gage_results$less %>% as.data.frame() %>% filter(q.val < pval_threshold)
  
  # Combine results
  significant_results <- bind_rows(
    greater_significant %>% mutate(Direction = "Upregulated"),
    less_significant %>% mutate(Direction = "Downregulated")
  ) %>% arrange(p.val)  # Sort by p-value for better readability
  
  # Save results in a single file
  results_file <- file.path(output_path, paste0("differentially_expressed_", control, "_vs_", sample, ".tsv"))
  write.table(significant_results, file = results_file, sep = "\t", quote = FALSE, row.names = TRUE)
  }
  else {
    results_file <- file.path(output_path, paste0("differentially_expressed_bi_", control, "_vs_", sample, ".tsv"))
    write.table(sig_results, file = results_file, sep = "\t", quote = FALSE, row.names = TRUE)
  }
  
  return(sig_results)
}

# Example usage:
# run_gage_analysis("path/to/feature_counts.tsv", "path/to/custom_pathways.gmt", "output_directory/gage_results")
