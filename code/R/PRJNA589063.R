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
library(dplyr)
library(readr)
library(gage)
library(tibble)
library(ggplot2)


featureCountsAbsolutePathGetmm = "C:/Users/Bisbii/PythonProjects/omics-integration/data/ngaditana/PRJNA589063/getmm_with_replicates.tsv"

getmm = read.delim(featureCountsAbsolutePathGetmm, header = TRUE, row.names = 1)



sample_names <- colnames(getmm)
Treat <- factor(substr(sample_names, 1, 1))
Time <- factor(substr(sample_names, 2, 2))

y <- DGEList(counts=getmm, group=Treat)



keep <- filterByExpr(y, min.count=1, min.total.count=1)
table(keep)

y <- y[keep, , keep.lib.sizes=FALSE]

#colors <- brewer.pal(nlevels(Treat), "Set1") # Choose a color palette
pdf("C:/Users/Bisbii/PythonProjects/omics-integration/results/ngaditana/PRJNA589063/mds_rnaseq_ng.pdf", width = 12)
plotMDS(y, col=colors[Treat], pch=c(15,16), cex = 1.5, cex.lab=1.3)
legend("topright", legend=levels(Treat), col=colors, pch=16, pt.bg=colors, cex=1.3, title="Condition")
legend("bottomright", legend=levels(Time), pch=c(16,15), col="black", cex=1.3, title="Time")
dev.off()


design <- model.matrix(~Time+Time:Treat)
row.names(design) = colnames(getmm)


logFC <- predFC(y,design,prior.count=1,dispersion=0.05)

cor(logFC)


Group <- factor(paste(Treat,Time,sep="."))

design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)

logFC <- predFC(y,design,prior.count=1,dispersion=0.05)

cor(logFC)


y <- estimateDisp(y, design, robust = TRUE)

y$common.dispersion

fit <- glmQLFit(y, design, robust=TRUE)


con_nacl <- makeContrasts(
  + C.3 - N.3, 
  + C.3 - P.3,
  + C.5 - N.5, 
  + C.5 - P.5,
  levels=Group)

#contrast = makeContrasts(C.5 - N.5,levels=Group)

qlf <- glmQLFTest(fit, contrast=con_nacl)

topTags(qlf)

summary(decideTests(qlf))

significant_genes <- topTags(qlf, n = Inf)$table
significant_genes <- significant_genes[significant_genes$FDR <= 0.05, ]
#significant_genes <- significant_genes[significant_genes$logFC > 0, ]
significant_gene_ids <- rownames(significant_genes)

aggregated_counts <- sumTechReps(getmm[significant_gene_ids, ], ID = Group)

heatmap(as.matrix(aggregated_counts))

write.table(significant_genes, file = "C:/Users/Bisbii/PythonProjects/omics-integration/results/ngaditana/PRJNA589063/deg/degs.tsv", 
            sep = "\t", quote = FALSE)


contrast_c5_n5 = makeContrasts(C.5 - P.5,levels=Group)


qlf <- glmQLFTest(fit, contrast=contrast_c5_n5)

topTags(qlf)

summary(decideTests(qlf))

significant_genes <- topTags(qlf, n = Inf)$table
significant_genes <- significant_genes[significant_genes$FDR <= 0.05, ]
significant_gene_ids <- rownames(significant_genes)

aggregated_counts <- sumTechReps(getmm[significant_gene_ids, ], ID = Group)

heatmap(as.matrix(aggregated_counts))

write.table(significant_genes, file = "C:/Users/Bisbii/PythonProjects/omics-integration/results/ngaditana/PRJNA589063/deg/degs_c5_p5.tsv", 
            sep = "\t", quote = FALSE)



getmm = read.delim(featureCountsAbsolutePathGetmm, header = TRUE, row.names = 1)

rownames(getmm) <- gsub("-RA_1$", "_RA", rownames(getmm))

rownames(getmm) <- gsub("-RB_1$", "_RB", rownames(getmm))

cn=colnames(getmm)

custom_pathways <- read.gmt("C:/Users/Bisbii/PythonProjects/omics-integration/results/ngaditana/PRJNA589063/deg/custom_pathways_KOG.gmt")


gsets_list <- custom_pathways %>%
  group_by(term) %>%
  summarise(genes = list(gene)) %>%
  deframe()  # Converts to a named list


control_3 = grep('C3',cn, ignore.case =T)
nitrogen_3=grep('N3',cn, ignore.case =T)


gse16873.kegg.p <- gage(getmm, gsets = gsets_list, ref = control_3, samp = nitrogen_3, compare = "unpaired", same.dir=T, rank.test = T) #


gse16873.kegg.sig<-sigGeneSet(gse16873.kegg.p, outname="gse16873.kegg")



source("C:/Users/Bisbii/Documents/R/gage.R")
source("C:/Users/Bisbii/Documents/R/gage_plot.R")

featureCountsAbsolutePathGetmm = "C:/Users/Bisbii/PythonProjects/omics-integration/data/ngaditana/PRJNA589063/getmm_with_replicates.tsv"
getmm = read.delim(featureCountsAbsolutePathGetmm, header = TRUE, row.names = 1)

rownames(getmm) <- gsub("-RA_1$", "_RA", rownames(getmm))
rownames(getmm) <- gsub("-RB_1$", "_RB", rownames(getmm))

gmt = "C:/Users/Bisbii/PythonProjects/omics-integration/results/ngaditana/PRJNA589063/deg/custom_pathways_KO.gmt"
output = "C:/Users/Bisbii/PythonProjects/omics-integration/results/ngaditana/PRJNA589063/deg/gage_ko"
res = run_gage_analysis(getmm, gmt, output, "C3", "N3")


plot_gage_results(res, "C:/Users/Bisbii/PythonProjects/omics-integration/results/ngaditana/PRJNA589063/ng_n3_deg.pdf",
                  title = "N3 vs C3", qval_cutoff=0.01)


res = run_gage_analysis(getmm, gmt, output, "C5", "N5")

plot_gage_results(res, "C:/Users/Bisbii/PythonProjects/omics-integration/results/ngaditana/PRJNA589063/ng_n5_deg.pdf",
                  title = "N5 vs C5", qval_cutoff=0.01)

res = run_gage_analysis(getmm, gmt, output, "C3", "P3")

plot_gage_results(res, "C:/Users/Bisbii/PythonProjects/omics-integration/results/ngaditana/PRJNA589063/ng_p3_deg.pdf",
                  title = "P3 vs C3",, qval_cutoff=0.01)

res = run_gage_analysis(getmm, gmt, output, "C5", "P5")

plot_gage_results(res, "C:/Users/Bisbii/PythonProjects/omics-integration/results/ngaditana/PRJNA589063/ng_p5_deg.pdf",
                  title = "P5 vs C5", qval_cutoff=0.01)

res = run_gage_analysis(getmm, gmt, output, "C3", "N3", bidirectional=F)
res = run_gage_analysis(getmm, gmt, output, "C5", "N5", bidirectional=F)
res = run_gage_analysis(getmm, gmt, output, "C3", "P3", bidirectional=F)
res = run_gage_analysis(getmm, gmt, output, "C5", "P5", bidirectional=F)


