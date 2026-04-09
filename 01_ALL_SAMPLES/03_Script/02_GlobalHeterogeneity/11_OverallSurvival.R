library(tidyverse)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(broom)
library(forcats)

set.seed(1234)

# -- bulkRNA deconvolution using BayesPrism -- 

#github:https://github.com/Danko-Lab/BayesPrism 
#find tutorial here: https://github.com/Danko-Lab/BayesPrism/blob/main/tutorial_deconvolution.html

#single cell reference dataset
# FOR NSCLC - Kim et al. 2020 (GSE131907) - https://doi.org/10.1038/s41467-020-16164-1 
# FOR Breast Cancer - Bassez et al. 2021 - (http://biokey.lambrechtslab.org/) 


# -- Overall Survival analysis --

# 1. Load NK summarized expression objects derived through BayesPrism and the TCGA

se_list <- list(
  LUAD = readRDS("path/to/se_LUAD_norm.rds"),
  LUSC = readRDS("path/to/se_LUSC_norm.rds"),
  BRCA = readRDS("path/to/se_BRCA_norm.rds")
)

# 2. load signatures to test

NKsig <- readLines("path/to/13_NK_sig.txt")

taNK_sig <- read.csv("path/to/taNK_markers.csv")
rownames(taNK_sig) <- taNK_sig$X
nontaNK <- rownames(taNK_sig[order(taNK_sig$avg_log2FC), ])[1:10]
taNK <- rownames(taNK_sig[order(-taNK_sig$avg_log2FC), ])[1:10]

sigs_genes <- list(NKsig = NKsig, taNK = taNK, nontaNK = nontaNK)


# 3. Survival Helper Function

run_signature_survival <- function(se_obj,
                                   genes,
                                   tumor_code,
                                   groups = c("low","high")) {
  
  clin <- as_tibble(colData(se_obj)) |>
    mutate(
      status = vital_status,
      days = if_else(status == 1,
                     days_to_death,
                     days_to_last_followup)
    )
  
  genes_use <- intersect(genes, rownames(se_obj))
  
  if(length(genes_use) < 3){
    stop("Too few genes from signature found in expression matrix")
  }
  
  expr <- assay(se_obj)[genes_use, , drop = FALSE]
  
  expr_scaled <- t(scale(t(expr)))
  
  score <- colMeans(expr_scaled, na.rm = TRUE)
  
  # avoid duplicated quantiles
  score <- jitter(score, factor = 1e-6)
  
  q <- unique(quantile(score, probs = seq(0,1,0.25)))
  
  group <- cut(score, q, include.lowest = TRUE)
  levels(group) <- c("low","low_mid","mid","high")
  
  surv_df <- clin |>
    mutate(signature_score = score,
           group = group) |>
    filter(group %in% groups) |>
    drop_na(days, status, group)
  
  fit <- survfit(Surv(days, status) ~ group, data = surv_df)
  surv_df$group <- factor(surv_df$group, levels = c("low","high"))
  cox <- coxph(Surv(days, status) ~ group, data = surv_df)
  
  return(list(
    fit = fit,
    cox = cox,
    data = surv_df
  ))
}

# 4. Run ALL Signatures × ALL Tumors

results <- list()

for(tumor in names(se_list)){
  
  message("Running tumor: ", tumor)
  
  results[[tumor]] <- list()
  
  for(sig in names(sigs_genes)){
    
    message("   Signature: ", sig)
    
    results[[tumor]][[sig]] <-
      run_signature_survival(
        se_obj = se_list[[tumor]],
        genes = sigs_genes[[sig]],
        tumor_code = tumor
      )
  }
}

results[["LUAD"]][["NKsig"]]
results[["LUAD"]][["taNK"]]
results[["BRCA"]][["nontaNK"]]


# 5. plot Kaplan Maier - Fig. 4F
dir.create("KM_plots", showWarnings = FALSE)

for(tumor in names(results)){
  for(sig in names(results[[tumor]])){
    
    fit <- results[[tumor]][[sig]]$fit
    df  <- results[[tumor]][[sig]]$data
    
    p <- ggsurvplot(
      fit,
      data = df,
      pval = TRUE,
      conf.int = FALSE,
      palette = c("#009999", "#FF6666"),
      risk.table = FALSE,
      title = paste(tumor, sig)
    )
    
    p2 <- p$plot +
      theme(
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
      )
    
    ggsave(
      filename = paste0("KM_plots/KM_", tumor, "_", sig, ".pdf"),
      plot = p2,
      width = 4,
      height = 4)
  }
}

#6. Extract HR

extract_forest_table <- function(results, signature){
  
  res <- list()
  
  for(tumor in names(results)){
    
    cox <- results[[tumor]][[signature]]$cox
    
    tmp <- tidy(cox) |>
      mutate(
        tumor = tumor,
        HR = exp(estimate),
        lower_ci = exp(estimate - 1.96 * std.error),
        upper_ci = exp(estimate + 1.96 * std.error)
      )
    
    res[[tumor]] <- tmp
  }
  
  bind_rows(res)
}

#7. Forest Plot Function - Fig. S8A

plot_forest <- function(table,
                        outfile){
  
  p <- ggplot(table,
              aes(x = fct_reorder(tumor, HR),
                  y = HR)) +
    geom_hline(yintercept = 1,
               linetype = "dashed",
               linewidth = 0.3) +
    geom_point(size = 3,
               color = "#FF6666") +
    geom_errorbar(aes(ymin = lower_ci,
                      ymax = upper_ci),
                  width = 0.2,
                  color = "#FF6666") +
    scale_y_log10() +
    coord_flip() +
    theme_minimal() +
    labs(y = "Hazard Ratio (95% CI)",
         x = "") +
    theme(
      plot.title = element_blank(),
      #axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size =11, colour = "black", face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, colour = "black"),
      axis.title.x = element_text(size = 10, face = "bold"),
      panel.grid =element_blank(),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.3),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line()
    )
  
  ggsave(outfile, p, width = 3, height = 3)
  
  return(p)
}


#8. Forrest Plots - all three sigs

forest_tables <- list()

for(sig in names(sigs_genes)){
  
  message("Forest plot for: ", sig)
  
  tab <- extract_forest_table(results, sig)
  
  write.csv(tab,
            paste0("HR_", sig, ".csv"),
            row.names = FALSE)
  
  plot_forest(tab,
              paste0("Forest_", sig, ".pdf"))
  
  forest_tables[[sig]] <- tab
}

#save
save(results,
     forest_tables,
     file = "NK_signature_survival_results.RData")

