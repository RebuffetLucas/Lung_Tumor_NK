library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(hdf5r)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

fragpath <- list()
counts <- list()
Seurats <- list()

#create Seurat Objects

for (i in 1:11){
  fragpath[[i]]<- paste("path/to/CR_outs_ARC_2.0.2/CSS",i,"_count/outs/atac_fragments.tsv.gz", sep="")
  names(fragpath)[[i]] <- paste("CSS", i, sep = "")
  counts[[i]] <- Read10X_h5(paste("path/to/CR_outs_ARC_2.0.2/CSS",i,"_count/outs/filtered_feature_bc_matrix.h5", sep = ""))
  names(counts)[[i]] <- paste("CSS", i, sep = "")
  Seurats[[i]] <- CreateSeuratObject(
    counts = counts[[i]]$"Gene Expression",
    assay = "RNA")
  Seurats[[i]][["ATAC"]] <- CreateChromatinAssay(
    counts = counts[[i]]$Peaks,
    sep = c(":", "-"),
    fragments = fragpath[[i]],
    annotation = annotation)
  
}

saveRDS(fragpath, "path/to/fragpath.rds")

## --PEAK CALLING (MACS2) --

for (i in seq_along(Seurats)){
  DefaultAssay(Seurats[[i]]) <- "ATAC"
  peaks <- CallPeaks(Seurats[[i]], macs2.path = "path/to/MACS2/2.2.7.1-foss-2021a/bin/macs2")
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  macs2_counts <- FeatureMatrix(
    fragments = Fragments(Seurats[[i]]),
    features = peaks,
    cells = colnames(Seurats[[i]])
  )
  Seurats[[i]][["peaks"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = fragpath[[i]],
    annotation = annotation
  )
}

#Obtain common peak set

combined.peaks <- reduce(x = do.call(c, lapply(Seurats, function(x) x@assays$peaks@ranges)))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10000 & peakwidths > 20]

new.counts <- list()

# Requantify unified set of peaks via FeatureMatrix()

for (i in seq_along(Seurats)){
  new.counts[[i]] <- FeatureMatrix(
    fragments = Seurats[[i]]@assays$peaks@fragments,
    features = combined.peaks,
    cells = rownames(Seurats[[i]]@meta.data))
  gc()
  Seurats[[i]][["peaks"]] <- CreateChromatinAssay(
    counts = new.counts[[i]],
    fragments = fragpath[[i]],
    annotation = annotation
  )
}

saveRDS(Seurats, "path/to/Seurats_list.RDS")





