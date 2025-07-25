library(Seurat)
library(Signac)
library(ggpubr)

Seurats <- readRDS("path/to/Seurats_list.RDS")
fragpath <- readRDS("path/to/fragpath.rds")

#----ATAC QC----

#calculate FRiP, blacklist regions, nucleosome signal, TSS per sample

# Load blacklist regions
blacklist_hg38 <- Signac::Blacklist(genome = "hg38")

for (i in seq_along(Seurats)) {
  DefaultAssay(Seurats[[i]]) <- "peaks"
  total_fragments <- CountFragments(fragpath[[i]])
  Seurats[[i]]$fragments <- total_fragments[match(colnames(Seurats[[i]]), total_fragments$CB), "frequency_count"]
  Seurats[[i]] <- FRiP(
    object = Seurats[[i]],
    assay = "peaks",
    total.fragments = "fragments"
  )
  DefaultAssay(Seurats[[i]]) <- "ATAC"
  Seurats[[i]]$blacklist_fraction <- FractionCountsInRegion(
    object = Seurats[[i]],
    assay = "ATAC",
    regions = blacklist_hg38
  )
  Seurats[[i]] <- NucleosomeSignal(Seurats[[i]])
  Seurats[[i]] <- TSSEnrichment(Seurats[[i]], fast = FALSE)
}

# QC Plotting
plt1 <- lapply(seq_along(Seurats), function(i) {
  Seurats[[i]]$nucleosome_group <- ifelse(Seurats[[i]]$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
  FragmentHistogram(object = Seurats[[i]], group.by = 'nucleosome_group') +
    ggtitle(Seurats[[i]]@meta.data$sample[1])
})

plt2 <- lapply(seq_along(Seurats), function(i) {
  Seurats[[i]]$high.tss <- ifelse(Seurats[[i]]$TSS.enrichment > 2.5, 'High', 'Low')
  TSSPlot(Seurats[[i]], group.by = 'high.tss') + NoLegend() +
    ggtitle(Seurats[[i]]@meta.data$sample[1])
})

pdf(file = "path/to/ATAC_Nucleosome.pdf", width = 20, height = 12)
print(ggarrange(plotlist = plt1, nrow = 3, ncol = 4))
dev.off()

pdf(file = "path/to/ATAC_TSS.pdf", width = 20, height = 12)
print(ggarrange(plotlist = plt2, nrow = 3, ncol = 4))
dev.off()

# Merge for RNA QC
sample <- c("CSS1", "CSS4", "CSS7", "CSS10",
            "CSS13", "CSS16", "CSS19", "CSS21",
            "CSS23", "CSS25", "CSS27")

merged <- merge(Seurats[[1]], Seurats[2:length(Seurats)],
                add.cell.ids = sample,
                merge.data = TRUE)

#----RNA QC----

DefaultAssay(merged) <- "RNA"

#save percentages for ribo, mito + stress genes
merged[["percent.rb"]] <- PercentageFeatureSet(merged, pattern = "^RP[S|L]")
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged[["percent.HSP"]] <- PercentageFeatureSet(merged, pattern ="^HSP")

#plotting
QC.features <- c("nFeature_RNA", "nCount_RNA","percent.mt",
                 "percent.rb", "percent.HSP")

pdf(file = "path/to/RNA_QC.pdf", width = 15, height = 12)
VlnPlot(merged, group.by = "sample", features = QC.features, pt.size = 0.1, ncol = 2) +
  NoLegend()
dev.off()

#CellCycle
merged <- CellCycleScoring(object = merged, g2m.features = cc.genes$g2m.genes,
                              s.features = cc.genes$s.genes)

pdf(file = "path/to/CellCycle.pdf", width = 6, height = 3)
VlnPlot(merged, features = c("S.Score", "G2M.Score"), group.by = "sample",
        ncol = 2, pt.size = 0.1)
dev.off()


#sample gender
gene_table <- read.csv("path/to/genes.table.csv")
#gene table downloaded from "https://github.com/NBISweden/workshop-scRNAseq/blob/master/labs/misc/genes.table.csv"
#head(gene_table)

chrY.gene = gene_table$external_gene_name[gene_table$chromosome_name == "Y"]
merged$pct_chrY = colSums(merged@assays$RNA@counts[chrY.gene, ])/colSums(merged@assays$RNA@counts)
FeatureScatter(merged, feature1 = "XIST", feature2 = "pct_chrY", group.by = "sample")

pdf(file = "path/to/gender.pdf", width = 8)
VlnPlot(merged, features = c("XIST", "pct_chrY"))
dev.off()

#subset dataset based on selected thresholds

merged_subset <- subset(
  x = merged,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nFeature_RNA > 200 &
    nFeature_ATAC > 200 &
    percent.mt < 20 &
    percent.rb > 0.05 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    blacklist_fraction < 0.05 &
    FRiP > 0.25
)


saveRDS(merged, "path/to/merged_subset.RDS")










