library(BayesPrism)
library(Seurat)
library(SummarizedExperiment)
library(TCGAutils)
library(curatedTCGAData)
library(tidyverse)
library(survival)
library(ggpubr)
library(survminer)
library(multcomp)
library(broom)

# -- bulkRNA deconvolution using BayesPrism -- 

#github:https://github.com/Danko-Lab/BayesPrism 
#find tutorial here: https://github.com/Danko-Lab/BayesPrism/blob/main/tutorial_deconvolution.html

#BayesPrism uses a reference single cell dataset with cell type information
#it can estimate cell type fraction and cell type specific gene expression Z in the bulk

#load single cell refernce dataset - here lung adenocarcinoma Kim et al. 2020 (GSE131907) - https://doi.org/10.1038/s41467-020-16164-1 
sc_Lung <- readRDS(gzfile("path/to/GSE131907_LUAD_raw_UMI_matrix.rds"))
metadata <- read.delim("path/to/GSE131907_LUAD_annotation.txt")

#add rownames to metadata to fit colnames in matrix
rownames(metadata) <- metadata$Index
intersect(rownames(metadata), colnames(sc_Lung)) #check

#get just tumor tissue
tumor_barcodes <- metadata$Index[metadata$Sample_Origin == "tLung"]
sc_tumor <- sc_Lung[,colnames(sc_Lung) %in% tumor_barcodes] #subset the matrix for barcodes only from tumor tissue

#subset metadata
met <- metadata[, c("Index","Cell_type", "Sample_Origin", "Sample")]
met <- subset(met, Sample_Origin == "tLung")
met$cell_state <- paste(met$Cell_type, met$Sample, sep = "_")

#prepare cell state and cell type labels - see vignette
cell.state.labels_LUAD <- met$cell_state
cell.type.labels_LUAD <- met$Cell_type

#convert into dense matrix + transpose, rownames should be cell IDs, colnames genes names
sc_tumor <- as.matrix(sc_tumor)
sc_tumor <- t(sc_tumor)

sort(table(cell.type.labels_LUAD)) #have a look

table(cbind.data.frame(cell.state.labels_LUAD, cell.type.labels_LUAD))

#QC - single cell reference

#pariwise correlation matrix between cell states and cell types - sense of their quality - see vignette
plot.cor.phi(input = sc_tumor, input.labels = cell.type.labels_LUAD, cexRow = 0.2, cexCol = 0.2)

#filter outlier genes --> mean expression across all cell types + cell type specificity scores
sc.stat <- plot.scRNA.outlier(sc_tumor, cell.type.labels = cell.type.labels_LUAD, species = "hs", return.raw = T)
head(sc.stat)

#remove genes such as ribosomal/mitochondrial genes/MALAT1 + genes from X- & Y- chromosome
sc_Lung_tumor_filtered <- cleanup.genes(input = sc_tumor, input.type = "count.matrix", species = "hs", 
                                        gene.group = c("Rb",
                                                       "Mrp","other_Rb", "chrM",
                                                       "MALAT1", "chrX","chrY"),
                                        exp.cells = 5)

save.image("path/to/Kim_LUAD_sc_tumor.RData") #save image

#now - load the raw TCGA bulk datatset for LUAD and LUSC

only_Lung <- readRDS("path/to/Lung_rawcounts.RDS")

#only primary tumor
primary_tumors_lung <- TCGAprimaryTumors(only_Lung)

#take RAW counts
se_LUAD <- getWithColData(primary_tumors_lung, "LUAD_RNASeqGene-20160128")
se_LUSC <- getWithColData(primary_tumors_lung, "LUSC_RNASeqGene-20160128")

#get matrix
mt.LUAD <- assays(se_LUAD)[[1]]
mt.LUSC <- assays(se_LUSC)[[1]]

#row names need to be sample IDs, colnames need to be genes - check
head(colnames(mt.LUAD))
head(rownames(mt.LUAD))

#transpose if needed
mt.LUAD <- t(mt.LUAD)
mt.LUSC <- t(mt.LUSC)


#QC - bulk

#visualize outlier genes in bulkRNA
LUAD.stat <- plot.bulk.outlier(
  bulk.input=mt.LUAD,
  sc.input=sc_Lung_tumor_filtered, 
  cell.type.labels=cell.type.labels_LUAD,
  species="hs", 
  return.raw=TRUE)

LUSC.stat <- plot.bulk.outlier(
  bulk.input=mt.LUSC,
  sc.input=sc_Lung_tumor_filtered, 
  cell.type.labels=cell.type.labels_LUAD,
  species="hs", 
  return.raw=TRUE)

#check the concordance of gene expression for different types of genes
plot.bulk.vs.sc(sc.input = sc_Lung_tumor_filtered,
                bulk.input = mt.LUAD)

#subset protein coding genes
sc_lung.filtered.pc <- select.gene.type(sc_Lung_tumor_filtered,
                                        gene.type = "protein_coding")

#Select marker genes for each cell type
diff.exp.stat <- get.exp.stat(sc.dat=sc_Lung_tumor_filtered[,colSums(sc_Lung_tumor_filtered>0)>3],
                              cell.type.labels=cell.type.labels_LUAD,
                              cell.state.labels=cell.state.labels_LUAD,
                              pseudo.count=0.1,
                              cell.count.cutoff=50, 
                              n.cores=4)

# colnames(diff.exp.stat$`NK cells`)
# head(diff.exp.stat$`NK cells`)

#subset count matrix over the signature genes
sc_lung.filtered.pc_sig <- select.marker(sc.dat=sc_lung.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.05,
                                         lfc.min=0.1)

#Contruct a prism object
myPrism_LUAD <- new.prism(
  reference=sc_lung.filtered.pc,
  mixture=mt.LUAD,
  input.type="count.matrix",
  cell.type.labels = cell.type.labels_LUAD,
  cell.state.labels = cell.state.labels_LUAD,
  key="Epithelial cells", #malignant cell type/tumor
  outlier.cut=0.01,
  outlier.fraction=0.1)

myPrism_LUSC <- new.prism(
  reference=sc_lung.filtered.pc,
  mixture=mt.LUSC,
  input.type="count.matrix",
  cell.type.labels = cell.type.labels_LUAD,
  cell.state.labels = cell.state.labels_LUAD,
  key="Epithelial cells", #malignant cell type/tumor
  outlier.cut=0.01,
  outlier.fraction=0.1)


#RUN PRISM

bp.res_LUSC  <- run.prism(prism = myPrism_LUSC, n.cores=50)
bp.res_LUAD <- run.prism(prism = myPrism_LUAD, n.cores=50)

# extract posterior mean of cell type fraction theta
theta_LUAD <- get.fraction (bp=bp.res_LUAD,
                            which.theta="final",
                            state.or.type="type")

theta_LUSC <- get.fraction (bp=bp.res_LUSC,
                            which.theta="final",
                            state.or.type="type")

#extract posterior mean of cell type-specific gene expression count matrix Z (NK cells)
NK_LUAD <- get.exp (bp=bp.res_LUAD,
                    state.or.type="type",
                    cell.name="NK cells")


NK_LUSC <- get.exp (bp=bp.res_LUSC,
                    state.or.type="type",
                    cell.name="NK cells")

#Cell Type specific gene expression Matrix


#normalization of bulkRNA + add to sumarized experiment
NK_LUAD <- vst(round(t(NK_LUAD)))
NK_LUSC <-  vst(round(t(NK_LUSC)))

x <- colData(se_LUAD)

rownames(colData(se_LUAD)) # check

se_LUAD_norm<- SummarizedExperiment(assays=list(counts=NK_LUAD),
                                    colData=colData(se_LUAD))

se_LUSC_norm <- SummarizedExperiment(assays=list(counts=NK_LUSC),
                                     colData=colData(se_LUSC))


save.image("path/to/LUAD_sc_calc.RData")

# -- Overall Survival analysis --

tNK <- read.csv("/path/to/taNK_markers.csv")
tNK <- tNK[order(tNK$avg_log2FC, decreasing = T),]
top10 <- tNK$X[1:10] #top10 genes

my_colors <- setNames(c("#009999","#CCC","#666","#FF6666"), c("low", "low_mid","mid","high"))
group <- c("low", "high")

se_NK <- list(LUAD=se_LUAD_norm, LUSC = se_LUSC_norm)

surv_data <- list()
gene <- list()
E <- list()
e_avg <- list()
e_group <- list()
surv_data <- list()
fit <- list()
stats <- list()
plt <- list()

for(i in 1:length(se_NK)){
  surv_data[[i]] <- colData(se_NK[[i]])[,getClinicalNames(names(se_NK)[[i]])] |>
    as_tibble() |>
    mutate(status = vital_status,
           days=if_else(status==1, days_to_death, days_to_last_followup))
  names(surv_data)[[i]] <- paste(names(se_NK)[[i]], "_surv", sep = "")
  gene[[i]] <- top10[top10 %in% rownames(se_NK[[i]])]
  names(gene)[[i]] <- paste(names(se_NK)[[i]], "_genes", sep = "")
  E[[i]] <- assay(se_NK[[i]])[gene[[i]],,drop=FALSE]
  E[[i]] <- t(scale(t(E[[i]]), center=TRUE, scale=TRUE))
  names(E)[[i]] <- paste(names(se_NK)[[i]], "_E", sep = "")
  e_avg[[i]] <-  colMeans(E[[i]], na.rm = TRUE)
  names(e_avg)[[i]] <- paste(names(se_NK)[[i]], "_eavg", sep = "")
  e_group[[i]] <- cut(e_avg[[i]], 
                      quantile(e_avg[[i]], seq(0,1, 0.25)), 
                      include.lowest = TRUE)
  names(e_group)[[i]] <- paste(names(se_NK)[[i]], "_egroup", sep = "")
  levels(e_group[[i]]) <- names(my_colors)
  surv_data[[i]] <- surv_data[[i]] |>
    mutate(gene_expression = e_avg[[i]],
           gene_group = e_group[[i]])
  surv_data[[i]] <- surv_data[[i]] |>
    dplyr::filter(gene_group %in% group) |>
    mutate(gene_group = droplevels(gene_group))
  surv_data[[i]] <- surv_data[[i]][complete.cases(surv_data[[i]]$status, 
                                                  surv_data[[i]]$days, 
                                                  surv_data[[i]]$gene_group),]
  fit[[i]] <- survfit(Surv(days, status) ~ gene_group, data = surv_data[[i]])
}

stats <- list()
for(i in 1:length(se_NK)){
  stats[[i]] <- coxph(Surv(days, status) ~ gene_group, data = surv_data[[i]])}


names(fit) <- names(se_NK)
names(stats) <- names(se_NK)

p <- ggsurvplot_list(fit, data = surv_data, 
                     size = 1.5,                
                     palette = c("#009999", "#FF6666"),
                     conf.int = F,         
                     pval = TRUE,           
                     pval.method = TRUE,
                     risk.table = F,
                     title = names(se_NK)) 

pl <- list()

for(i in 1:length(p)){
  pl[[i]] <- p[[2]]$plot +  theme(axis.title.x = element_text(size = 12),
                                  axis.title.y = element_text(size = 15, face = "bold"),
                                  plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
}


#plot KM curves - Fig. 3G + Fig. S6B
pdf("path/to/KM_curve_LUAD.pdf", 
    width =3.8, height =4)
p[[1]]$plot +  theme(axis.title.x = element_text(size = 12),
                     axis.title.y = element_text(size = 14, face = "bold"),
                     plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
dev.off()

pdf("path/to/KM_curve_LUSC.pdf", 
    width =3.8, height =4)
p[[2]]$plot +  theme(axis.title.x = element_text(size = 12),
                     axis.title.y = element_text(size = 14, face = "bold"),
                     plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
dev.off()


# Forrest Plot - Fig. S6A
names(stats) <- names(se_NK)

cm.glht <- list()
pval_tab <- list()

for (i in 1:length(stats)){
  cm.glht[[i]] <- glht(stats[[i]], mcp(gene_group = "high - low = 0"))
  names(cm.glht)[[i]] <- names(stats)[[i]]
  pval_tab[[i]] <- tidy(cm.glht[[i]]) |>
    dplyr::select(-c(term, null.value))
  names(pval_tab)[[i]] <- names(stats)[[i]]
  pval_tab[[i]][["tumor_type"]] <- names(pval_tab)[[i]]
}


for (i in 1:length(surv_data)){
  surv_data[[i]][["tumor_type"]] <- names(se_NK)[[i]]
}


names(surv_data) <- names(se_NK)
surv_data_sub <- surv_data[c(names(pval_tab))]

table_p <- data.table::rbindlist(pval_tab)
table_surv <- data.table::rbindlist(surv_data_sub, fill = T)

table_surv$tumor_type <- as.factor(table_surv$tumor_type)
table_p$tumor_type <- as.factor(table_p$tumor_type)

#calculate HR
table_p <- table_p %>%
  mutate(
    HR = exp(estimate),  # Transforming estimate to hazard ratio
    lower_ci = exp(estimate - 1.96 * std.error),
    upper_ci = exp(estimate + 1.96 * std.error)
  )

#FORREST PLOT
pdf("path/to/ForrestPlot.pdf",width = 2.5, height = 2)
ggplot(table_p, aes(x = fct_reorder(tumor_type, HR), y = HR)) + 
  geom_hline(yintercept = 1, linetype = "dashed",linewidth = 0.3) +
  geom_point(size = 3, color = "#FF6666") +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, color = "#FF6666") +
  scale_y_log10() + 
  theme_minimal() +
  labs(
    y = "Hazard Ratio (95% CI)"
  )+
  coord_flip() +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size =11, colour = "black", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, colour = "black"),
    axis.title.x = element_text(size = 10, face = "bold"),
    panel.grid =element_blank(),
    axis.line.y.left = element_line(colour = "black", linewidth = 0.3),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.3),
    axis.ticks = element_line())
dev.off()

save.image("path/to/Survival_analysis.RData")
