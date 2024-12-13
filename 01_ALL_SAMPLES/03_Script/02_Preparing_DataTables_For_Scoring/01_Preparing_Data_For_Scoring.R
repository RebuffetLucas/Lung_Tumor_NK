#Preparing data for scoring for Clara
#This code aims at extracting a series of reference signatures that will be used in the project
library(openxlsx)

FILTER_pvalue_THR = 0.05
NUMBER_TOP_SCORING= 20


#Malmberg
markers= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Malmberg.rds")
write.xlsx(markers, "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Lung_Cancer_Basel/Gene_SIGNATURES_For_Clara/Malmberg.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = FALSE)



#Meta-NK TOP20 / TOP40

Markers_Seurat= readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/DEG/All_Markers_CITEseq3clusters.rds") #Already ready for scoring

NUMBER_TOP_SCORING= 20

Markers_Seurat %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC)-> top_All

write.xlsx(top_All, "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Lung_Cancer_Basel/Gene_SIGNATURES_For_Clara/MetaNK_CITEseq_Top20.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = FALSE)



NUMBER_TOP_SCORING= 40

Markers_Seurat %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  #filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC)-> top_All

write.xlsx(top_All, "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Lung_Cancer_Basel/Gene_SIGNATURES_For_Clara/MetaNK_CITEseq_Top40.xlsx", sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = FALSE)


#Convert GG signatures from mouse to human
FILE_SIGNATURES_PATH="/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK2/FromJanine/01_Reference/02_MonitoredGenesList/Ref_GG_Mouse_resident_circu.xlsx"


read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}



#Monitored markers
MONITORED_Markers = read_excel_allsheets(FILE_SIGNATURES_PATH)


MONITORED_Markers$Residency_score -> Residency_genes
MONITORED_Markers$Circulating_score -> Circulating_genes


  #Def functions
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  # ensembl <- useMart(biomart = "ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  mouse_attributes <- c("external_gene_name","external_gene_source","hsapiens_homolog_associated_gene_name","hsapiens_homolog_perc_id",
                        "hsapiens_homolog_goc_score", "hsapiens_homolog_orthology_confidence")
  genesV2<-getLDS(attributes = "hgnc_symbol",
                  filters = "hgnc_symbol", values = x, mart = ensembl.human,
                  attributesL = mouse_attributes, martL = ensembl.mouse,uniqueRows = TRUE)
  mousex <- genesV2 %>%
    filter(Gene.name != "") %>%
    group_by(HGNC.symbol) %>%
    slice_max(Human.Gene.order.conservation.score) %>%
    filter(Human.orthology.confidence..0.low..1.high. == 1) %>%
    dplyr::select(HGNC.symbol,Gene.name) %>%
    dplyr::rename(human.gene.name = HGNC.symbol, mouse.gene.name = Gene.name)
  
  return(mousex)
}


convertMouseGeneList <- function(x){
  
  require("biomaRt")
  # ensembl <- useMart(biomart = "ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  human_attributes <- c("external_gene_name","external_gene_source","mmusculus_homolog_associated_gene_name","mmusculus_homolog_perc_id",
                        "mmusculus_homolog_goc_score", "mmusculus_homolog_orthology_confidence")
  genesV2<-getLDS(attributes = "mgi_symbol",
                  filters = "mgi_symbol", values = x, mart = ensembl.mouse,
                  attributesL = human_attributes, martL = ensembl.human,uniqueRows = TRUE)
  
  humanx <- genesV2 %>%
    filter(Gene.name != "") %>%
    group_by(MGI.symbol) %>%
    slice_max(Mouse.Gene.order.conservation.score) %>%
    filter(Mouse.orthology.confidence..0.low..1.high. == 1) %>%
    dplyr::select(MGI.symbol,Gene.name) %>%
    dplyr::rename(mouse.gene.name = MGI.symbol, human.gene.name = Gene.name)
  
  return(humanx)
}

Human_eq_residency= convertMouseGeneList(Residency_genes)

Human_eq_residency <- Human_eq_residency %>%
  left_join(Human_eq_residency, by = c("gene" = "mouse.gene.name"))

#write.csv(Human_eq_residency , "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Lung_Cancer_Basel/Gene_SIGNATURES_For_Clara/GG_Residency_MouseToHuman.csv", row.names=FALSE) #Save for Data sup

Human_eq_Circulating= convertMouseGeneList(Circulating_genes)


Human_eq_Circulating <- Human_eq_Circulating %>%
  left_join(Human_eq_Circulating, by = c("gene" = "mouse.gene.name"))

#write.csv(Human_eq_Circulating , "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Lung_Cancer_Basel/Gene_SIGNATURES_For_Clara/GG_Circulating_MouseToHuman.csv", row.names=FALSE) #Save for Data sup



