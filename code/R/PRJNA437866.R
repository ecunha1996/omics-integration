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


featureCountsAbsolutePathGetmm = "C:/Users/Bisbii/PythonProjects/omics-integration/data/dsalina/PRJNA437866/getmm_with_replicates.tsv"

getmm = read.delim(featureCountsAbsolutePathGetmm, header = TRUE, row.names = 1)


sample_names <- colnames(getmm)


condition <- factor(sapply(sample_names, function(x) substr(x, 1, nchar(x) - 2)))

dge <- DGEList(counts = getmm, group = condition)


keep <- filterByExpr(dge)
summary(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]


mylegend= c(levels(condition))

mylegend[mylegend == "control"] <- "Control"
mylegend[mylegend == "nacl"] <- "NaCl"
mylegend[mylegend == "sorb"] <- "Sorbitol"
mylegend[mylegend == "h2o2"] <- expression(H[2]*O[2])


colors <- brewer.pal(nlevels(condition), "Set1") # Choose a color palette
pdf("C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA437866/mds_rnaseq_ds_stress.pdf")
mds = plotMDS(dge, col=colors[condition], pch=c(16), dim.plot = c(1,2), cex = 2, cex.lab=1.3)
legend("bottomright", legend=mylegend, col=colors, pch=c(16), pt.bg=colors, cex=1.3, title="Condition")
dev.off()

pdf("C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA437866/mds_rnaseq_ds_stress_1v3.pdf")
mds = plotMDS(dge, col=colors[condition], pch=c(16), dim.plot = c(1,3), cex = 2, cex.lab=1.3)
legend("bottomleft", legend=mylegend, col=colors, pch=c(16), pt.bg=colors, cex=1.3, title="Condition")
dev.off()



library(edgeR)

# Calculate logFC using edgeR
design <- model.matrix(~ condition)  # Replace 'condition' with your factor variable
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)




design <- model.matrix(~ 0 + condition)

colnames(design) <- levels(condition)

dge <- estimateDisp(dge, design, robust = TRUE)

dge$common.dispersion

plotBCV(dge)

fit <- glmQLFit(dge, design, robust=TRUE)


con_nacl <- makeContrasts(
  + nacl - control, 
  + h2o2 - control,
  + sorb - control,
  levels=design)

qlf <- glmQLFTest(fit, contrast=con_nacl)


topTags(qlf)

summary(decideTests(qlf))

significant_genes <- topTags(qlf, n = Inf)$table
significant_genes <- significant_genes[significant_genes$PValue < 0.05, ]
significant_gene_ids <- rownames(significant_genes)

aggregated_counts <- sumTechReps(getmm[significant_gene_ids, ], ID = condition)

heatmap(as.matrix(aggregated_counts))

write.table(significant_genes, file = "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA437866/deg/degs.tsv", 
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
      valid_contrasts <- c("nacl - control","h2o2 - control", "sorb - control")
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
        
        # custom_pathways <- custom_pathways[!is.na(custom_pathways$gene) & !is.na(custom_pathways$term), ]
        
        
        num_matching_genes <- sum(names(gene_list_func) %in% custom_pathways$gene)
        
        
        message("Genes matching pathways: ", num_matching_genes)
        
        
        gsea_results <- GSEA(geneList = gene_list_func, 
                             TERM2GENE = custom_pathways, 
                             minGSSize = 1, 
                             pvalueCutoff = 0.05,
                             maxGSSize = 50000,
                             seed=TRUE)
        
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
group_levels <- c("h2o2", "nacl", "sorb",  "control")  # Adjust based on your dataset

# Define path to custom pathways
gmt_file <- "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA495151/deg/custom_pathways_KO.gmt"


output_dir <- "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA437866/deg/gsea_results_KO"
# Run GSEA for all contrasts
gsea_results_all <- run_gsea_for_all_conditions(fit, group_levels, gmt_file, output_dir, minGSSize=5, pvalueCutoff=0.05)


source("C:/Users/Bisbii/Documents/R/gage.R")
source("C:/Users/Bisbii/Documents/R/gage_plot.R")
featureCountsAbsolutePathGetmm = "C:/Users/Bisbii/PythonProjects/omics-integration/data/dsalina/PRJNA437866/getmm_with_replicates.tsv"
getmm = read.delim(featureCountsAbsolutePathGetmm, header = TRUE, row.names = 1)
rownames(getmm) <- gsub("_", ".", rownames(getmm))
gmt = "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA437866/deg/custom_pathways_KO.gmt"
output = "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA437866/deg/gage_ko"


res = run_gage_analysis(getmm, gmt, output, "control", "h2o2")

plot_gage_results(res, 
                  "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA437866/ds_stress_h2o2_deg.pdf", 
                  title = "H2O2 vs Control", qval_cutoff=0.01)

res = run_gage_analysis(getmm, gmt, output, "control", "nacl")

plot_gage_results(res, 
                  "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA437866/ds_stress_nacl_deg.pdf", 
                  title = "NaCl vs Control", qval_cutoff=0.01)

res = run_gage_analysis(getmm, gmt, output, "control", "sorb")

plot_gage_results(res, 
                  "C:/Users/Bisbii/PythonProjects/omics-integration/results/dsalina/PRJNA437866/ds_stress_sorb_deg.pdf", 
                  title = "Sorbitol vs Control", qval_cutoff=0.01)

res = run_gage_analysis(getmm, gmt, output, "control", "h2o2", bidirectional=F)
res = run_gage_analysis(getmm, gmt, output, "control", "nacl", bidirectional=F)
res = run_gage_analysis(getmm, gmt, output, "control", "sorb", bidirectional=F)
