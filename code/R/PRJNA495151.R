library(pheatmap)
library("limma")
library("edgeR")
library("statmod")
library("genefilter")
library("Biobase")
library("readxl")
library("gplots")
library("TCC")
library(RColorBrewer)

featureCountsAbsolutePathGetmm = "C:/Users/Bisbii/PythonProjects/omics-integration/data/dsalina/PRJNA495151/getmm_with_replicates.tsv"

getmm = read.delim(featureCountsAbsolutePathGetmm, header = TRUE, row.names = 1)


sample_names <- colnames(getmm)

condition <- factor(sapply(sample_names, function(x) substr(x, 1, nchar(x) - 2)))

dge <- DGEList(counts = getmm, group = condition)


keep <- filterByExpr(dge, min.count=1, min.total.count=1)
summary(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]


colors <- brewer.pal(nlevels(condition), "Set1") # Choose a color palette
pdf("C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA495151/mds_rnaseq_ds_light.pdf", width = 12)
plotMDS(dge, col=colors[condition], pch=c(16), cex.lab=1.3, cex=1.5)
legend("topright", legend=levels(condition), col=colors, pch=16, pt.bg=colors, cex=1.5, title="Condition")
dev.off()


design <- model.matrix(~ 0 + condition)

colnames(design) <- levels(condition)

dge <- estimateDisp(dge, design, robust = TRUE)

dge$common.dispersion

plotBCV(dge)

fit <- glmQLFit(dge, design, robust=TRUE)


con_nacl <- makeContrasts(
  + ML - LL, 
  + HL - LL,
  levels=design)

qlf <- glmQLFTest(fit, contrast=con_nacl)


topTags(qlf)

summary(decideTests(qlf))

significant_genes <- topTags(qlf, n = Inf)$table
significant_genes <- significant_genes[significant_genes$PValue < 0.05, ]
significant_gene_ids <- rownames(significant_genes)

aggregated_counts <- sumTechReps(getmm[significant_gene_ids, ], ID = condition)

heatmap(as.matrix(aggregated_counts))

write.table(significant_genes, file = "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA495151/deg/degs.tsv", 
            sep = "\t", quote = FALSE)

run_gsea_for_all_conditions <- function(fit, group_levels, gmt_file=NULL, output_dir, minGSSize = 1, pvalueCutoff = 0.05, kog_file=NULL, ko_file=NULL) {
  library(edgeR)
  library(clusterProfiler)
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Load custom pathway database
  if (!is.null(gmt_file)){
    custom_pathways <- read.gmt(gmt_file)
  }
  
  # Store results for each contrast
  gsea_results_list <- list()
  
  # Loop through all pairwise comparisons
  for (i in 1:(length(group_levels) - 1)) {
    for (j in (i + 1):length(group_levels)) {
      valid_contrasts <- c("ML - LL", "HL - LL")
      contrast_name <- paste(group_levels[i], "-", group_levels[j], sep = "")
      message("Processing contrast: ", contrast_name)
      
      # Create contrast
      contrast_name = paste0(group_levels[i], " - ", group_levels[j])
      
      if (contrast_name %in% valid_contrasts) {
        contrast_matrix <- makeContrasts(contrasts = contrast_name ,levels=condition)
        
        # Perform differential expression analysis
        qlf <- glmQLFTest(fit, contrast = contrast_matrix)
        
        # Extract significant genes
        significant_genes <- topTags(qlf, n = Inf)$table
        
        rownames(significant_genes) <- paste0(rownames(significant_genes), ".1")
        
        if (nrow(significant_genes) == 0) {
          message("No significant genes found for ", contrast_name, ". Skipping.")
          next
        }
        
        gene_list_func <- significant_genes$logFC
        
        gene_list_func <- as.numeric(gene_list_func)
        
        
        gene_list_func <- sort(gene_list_func, decreasing = TRUE)
      
        names(gene_list_func) <- rownames(significant_genes)
        
        # Check how many genes match pathways
        
        custom_pathways <- custom_pathways[!is.na(custom_pathways$gene) & !is.na(custom_pathways$term), ]
        
        
        num_matching_genes <- sum(names(gene_list_func) %in% custom_pathways$gene)
        
        
        message("Genes matching pathways: ", num_matching_genes)
        
        
        gsea_results <- GSEA(geneList = gene_list_func, 
                             TERM2GENE = custom_pathways, 
                             minGSSize = minGSSize, 
                             pvalueCutoff = pvalueCutoff,
                             maxGSSize = 50000)
        
        # Save results
        gsea_results_list[[contrast_name]] <- gsea_results
        
        # Save to file
        output_file <- file.path(output_dir, paste0("gsea_results_", contrast_name, ".tsv"))
        write.table(gsea_results@result, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
        
        message("Finished processing: ", contrast_name)
      }
    }
  }
  
  return(gsea_results_list)
}

set.seed(42)

# Define conditions from your experiment
group_levels <- c("HL", "ML", "LL")  # Adjust based on your dataset

# Define path to custom pathways
gmt_file <- "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA495151/deg/custom_pathways_KO.gmt"


output_dir <- "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA495151/deg/gsea_results"
# Run GSEA for all contrasts
gsea_results_all <- run_gsea_for_all_conditions(fit, group_levels, gmt_file, output_dir, minGSSize=2, pvalueCutoff=0.05)


gsea_results_all$`HL - LL`@result

gsea_results_all$`ML - LL`@result


source("C:/Users/Bisbii/Documents/R/gage.R")
source("C:/Users/Bisbii/Documents/R/gage_plot.R")
featureCountsAbsolutePathGetmm = "C:/Users/Bisbii/PythonProjects/omics-integration/data/dsalina/PRJNA495151/getmm_with_replicates.tsv"
getmm = read.delim(featureCountsAbsolutePathGetmm, header = TRUE, row.names = 1)
rownames(getmm) <- gsub("_", ".", rownames(getmm))
gmt = "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA495151/deg/custom_pathways_KO.gmt"
output = "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA495151/deg/gage_ko"


res  = run_gage_analysis(getmm, gmt, output, "LL", "HL")

plot_gage_results(res, 
                  "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA495151/ds_light_hl_deg.pdf", 
                  title = "HL vs LL", qval_cutoff=0.01)

res  = run_gage_analysis(getmm, gmt, output, "LL", "ML")

plot_gage_results(res, 
                  "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA495151/ds_light_ml_deg.pdf", 
                  title = "ML vs LL", qval_cutoff=0.01)



pathways = read.gmt(gmt)

genes_in_ps = pathways$gene[pathways$term == "Citrate cycle (TCA cycle)"]

LL_samples <- grep("^LL_", colnames(getmm), value = TRUE)
HL_samples <- grep("^HL_", colnames(getmm), value = TRUE)

# Step 2: Compute mean expression for LL and HL
LL_mean <- rowMeans(getmm[, LL_samples])
HL_mean <- rowMeans(getmm[, HL_samples])

# Step 3: Compute log2 Fold Change (log2FC)
log2FC_vector <- log2(HL_mean / LL_mean)

# Name the vector by gene names
names(log2FC_vector) <- rownames(getmm)


upregulated_genes <- genes_in_ps[genes_in_ps %in% names(log2FC_vector[log2FC_vector > 0])]

downregulated_genes <- genes_in_ps[genes_in_ps %in% names(log2FC_vector[log2FC_vector < 0])]

res  = run_gage_analysis(getmm, gmt, output, "LL", "HL", bidirectional=F)
res  = run_gage_analysis(getmm, gmt, output, "LL", "ML", bidirectional=F)
