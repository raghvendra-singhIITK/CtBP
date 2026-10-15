# DESeq2 Analysis
# Differential gene expression and Targeted pathway heatmaps
rm(list = ls())
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(ashr)
library(showtext)
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
arrange <- dplyr::arrange
slice <- dplyr::slice
font_add("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_auto()
RESULT_FOLDER <- "Results"
OUT_DIR <- file.path(RESULT_FOLDER, "DESeq2_Pathway")
# DESeq2 thresholds
PADJ_THRESH <- 0.05
LFC_THRESH  <- 1.0
CELL_WIDTH  <- 0.42
FIG_HEIGHT  <- 3.4
FIG_DPI     <- 300

# Gene sets analyzed
# Aerobic Glycolysis
glycolysis_genes <- c("SLC2A1", "SLC2A2", "SLC2A3", "SLC2A4", "HK1", "HK2", "HK3", "GCK", "GPI", "PFKL", "PFKM", "PFKP", "PFKFB1", "PFKFB2", "PFKFB3", "PFKFB4",
  "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "PGK1", "PGK2", "PGAM1", "PGAM2", "ENO1", "ENO2", "ENO3", "PKM", "PKLR", "LDHA", "LDHB", "LDHC",
  "PDK1", "PDK2", "PDK3", "PDK4", "SLC16A1", "SLC16A3", "SLC16A7", "SLC16A8", "PHGDH", "PTBP1", "PTBP2", "MYC", "HIF1A")
# Hypoxia Response
hypoxia_genes <- c("HIF1A", "EPAS1", "ARNT", "ARNT2", "VHL", "EGLN1", "EGLN2", "EGLN3", "HIF1AN", "CUL2", "RBX1", "TCEB1", "TCEB2", "VEGFA", "VEGFB", "VEGFC",
  "ANGPT1", "ANGPT2", "SLC2A1", "LDHA", "PDK1", "PGK1", "ENO1", "ALDOA", "EPO", "TFRC", "TF", "CA9", "CA12", "SOD2", "GPX1", "TXNIP", "HMOX1",
  "BNIP3", "BNIP3L", "PINK1", "MTOR", "PRKAA1", "PRKAA2", "DDIT4", "LOX", "LOXL2", "MMP2", "MMP9", "P4HA1", "P4HA2", "CXCR4", "EGFR", "SERPINE1",
  "ADM", "NDRG1", "CITED2", "MXI1", "NOS2", "NOXA1", "RORA", "CHCHD4", "ATF4", "DDIT3", "XBP1", "CDKN1A")
# CtBP/EMT-invasion
ctbp_genes <- c("CTBP1", "CTBP2", "HDAC1", "HDAC2", "KDM1A", "RCOR1", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "CDH1", "CDKN2A", "PTEN", "BRCA1", "BRCA2",
  "MLH1", "TP53", "TP73", "PERP", "BAX", "BCL6", "RB1", "RBL1", "RBL2", "E2F1", "E2F2", "E2F3", "MDM2", "EZH2", "SUZ12", "BMI1", "RING1", "TRIM28",
  "ARID1A", "HIPK2", "RUNX1", "ZBTB16", "ZNF217", "NADK", "NNT", "VIM", "FN1", "CDH2", "LDHA", "LDHB", "GAPDH")

#Load data
config_path <- file.path(RESULT_FOLDER, "analysis_config.rds")
if (!file.exists(config_path))
  stop("analysis_config.rds not found at: ", config_path)
config <- readRDS(config_path)
tissues <- config$tissues
conditions <- config$conditions
meta <- config$metadata
cat("Cell lines:", paste(tissues,    collapse = ", "), "\n")
cat("Conditions:", paste(conditions, collapse = ", "), "\n")
cat("Samples   :", nrow(meta), "\n\n")
dir.create(OUT_DIR,                                 recursive = TRUE,  showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "Tier1_Full"),        recursive = FALSE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "Tier2_Significant"), recursive = FALSE, showWarnings = FALSE)
# Run DESeq2
all_results <- list()
for (tissue in tissues)
{
  cat(sprintf(">> %s\n", tissue))
  dds_path <- file.path(RESULT_FOLDER, tissue, paste0(tissue, "_dds.rds"))
  if (!file.exists(dds_path))
  {
    warning("DDS not found", tissue, ": ", dds_path)
    next
  }
  dds <- readRDS(dds_path)
  if (!"replicate" %in% colnames(colData(dds)))
  {
    rep_vals <- sub("^(Nx|Hx)(\\d+)$", "\\2", colnames(dds))
    colData(dds)$replicate <- factor(rep_vals)
    cat("replicate reconstructed from sample names\n")
  }
  colData(dds)$condition <- relevel(factor(colData(dds)$condition, levels = conditions), ref = conditions[1])
  design(dds) <- ~ replicate + condition
  dds <- DESeq(dds, test = "Wald", fitType = "parametric", quiet = TRUE)
  # ashr LFC shrinkage
  res <- lfcShrink(dds, contrast = c("condition", conditions[2], conditions[1]), type = "ashr", quiet = TRUE)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj) %>%
    mutate(cell_line = tissue,
           direction = case_when(!is.na(padj) & padj < PADJ_THRESH & log2FoldChange >  LFC_THRESH ~ "Up",
             !is.na(padj) & padj < PADJ_THRESH & log2FoldChange < -LFC_THRESH ~ "Down",TRUE ~ "NS"))
  all_results[[tissue]] <- res_df
  n_up <- sum(res_df$direction == "Up",   na.rm = TRUE)
  n_down <- sum(res_df$direction == "Down", na.rm = TRUE)
  cat(sprintf("   Up: %d  |  Down: %d  |  Total DEGs: %d\n\n", n_up, n_down, n_up + n_down))
  write.csv(res_df, file.path(OUT_DIR, paste0(tissue, "_DESeq2_Results.csv")), row.names = FALSE)
}
all_combined <- bind_rows(all_results)
write.csv(all_combined, file.path(OUT_DIR, "All_CellLines_DESeq2_Results.csv"), row.names = FALSE)
build_lfc_matrix <- function(all_results, gene_set, cell_lines)
{
  lfc_mat  <- matrix(NA_real_, nrow = length(cell_lines), ncol = length(gene_set), dimnames = list(cell_lines, gene_set))
  padj_mat <- matrix(1,        nrow = length(cell_lines), ncol = length(gene_set), dimnames = list(cell_lines, gene_set))
  for (cl in cell_lines)
  {
    df <- all_results[[cl]]
    if (is.null(df)) next
    df_sub <- df[df$gene %in% gene_set & !duplicated(df$gene), ]
    if (nrow(df_sub) == 0) next
    rownames(df_sub) <- df_sub$gene
    g <- intersect(gene_set, rownames(df_sub))
    lfc_mat[cl, g]  <- df_sub[g, "log2FoldChange"]
    padj_mat[cl, g] <- ifelse(is.na(df_sub[g, "padj"]), 1, df_sub[g, "padj"])
  }
  detected <- colSums(!is.na(lfc_mat)) > 0
  if (any(!detected))
    message("Not detected",
            paste(names(detected)[!detected], collapse = ", "))
  list(lfc = lfc_mat[, detected, drop = FALSE], padj = padj_mat[, detected, drop = FALSE])
}
compute_lfc_cap <- function(lfc_mat, padj_mat, padj_thresh = PADJ_THRESH, lfc_thresh = LFC_THRESH)
{
  sig_mask <- !is.na(padj_mat) & padj_mat < padj_thresh & abs(lfc_mat) >= lfc_thresh
  sig_lfcs <- abs(lfc_mat[sig_mask])
  if (length(sig_lfcs) == 0)
  {
    message("No significant genes")
    return(2)
  }
  cap <- max(ceiling(quantile(sig_lfcs, 0.95, na.rm = TRUE)), 2)
  cat(sprintf("   LFC cap (95th pct of sig genes): %d\n", cap))
  cap
}
filter_sig_genes <- function(lfc_mat, padj_mat,padj_thresh = PADJ_THRESH, lfc_thresh = LFC_THRESH)
{
  sig_mask <- !is.na(padj_mat) & padj_mat < padj_thresh & abs(lfc_mat) >= lfc_thresh
  keep     <- colSums(sig_mask, na.rm = TRUE) > 0
  if (sum(keep) == 0)
  {
    warning("No significant DEGs found")
    return(NULL)
  }
  cat(sprintf("   Significant genes (>= 1 cell line): %d / %d\n", sum(keep), ncol(lfc_mat)))
  for (g in colnames(lfc_mat)[keep])
  {
    sig_cls  <- rownames(lfc_mat)[sig_mask[, g]]
    lfc_vals <- round(lfc_mat[sig_cls, g], 2)
    cat(sprintf("     %-12s  %s  (log2FC: %s)\n", g, paste(sig_cls, collapse = ", "), paste(lfc_vals, collapse = ", ")))
  }
  list(lfc = lfc_mat[, keep, drop = FALSE], padj = padj_mat[, keep, drop = FALSE])
}
# Heatmaps
draw_pathway_heatmap <- function(lfc_mat, padj_mat, panel_title, panel_color, out_path, subtitle = NULL, lfc_cap,
                                 fig_height = FIG_HEIGHT, cell_width = CELL_WIDTH, dpi = FIG_DPI)
{
  n_genes <- ncol(lfc_mat)
  n_cl    <- nrow(lfc_mat)
  if (n_genes == 0)
  {
    warning("No genes", panel_title)
    return(invisible(NULL))
  }
  lfc_plot <- lfc_mat
  lfc_plot[lfc_plot >  lfc_cap] <-  lfc_cap
  lfc_plot[lfc_plot < -lfc_cap] <- -lfc_cap
  df <- as.data.frame(lfc_plot) %>%
    rownames_to_column("cell_line") %>%
    pivot_longer(-cell_line, names_to = "gene", values_to = "lfc") %>%
    left_join(as.data.frame(padj_mat) %>%
                rownames_to_column("cell_line") %>%
                pivot_longer(-cell_line, names_to = "gene", values_to = "padj_val"),
              by = c("cell_line", "gene")) %>%
    mutate(cell_line = factor(cell_line, levels = rev(rownames(lfc_mat))),
           gene = factor(gene, levels = colnames(lfc_mat)))
  clrs   <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF",
              "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
  breaks <- seq(-lfc_cap, lfc_cap, length.out = length(clrs))
  full_title <- if (!is.null(subtitle)) paste0(panel_title, "\n", subtitle) else panel_title
  p <- ggplot(df, aes(x = gene, y = cell_line)) +
    geom_tile(aes(fill = lfc), colour = "black", linewidth = 1, na.rm = TRUE) +
    geom_tile(data = filter(df, is.na(lfc)), fill = "#D0D0D0", colour = "black", linewidth = 0.5) +
    scale_fill_gradientn(colours  = clrs,
                         values   = scales::rescale(breaks, to = c(0, 1)),
                         limits   = c(-lfc_cap, lfc_cap),
                         oob      = scales::squish,
                         na.value = "#D0D0D0",
                         name     = expression(log[2]*FC),
                         guide    = guide_colorbar(barwidth = 1, barheight = 4.5,
                                                   title.position = "top", title.hjust = 0.5,
                                                   ticks.colour = "black", frame.colour = "black",
                                                   frame.linewidth = 0.4)) +
    scale_x_discrete(position = "bottom", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 13) +
    theme(axis.text.y      = element_text(family = "Arial", size = 30, face = "plain", colour = "black", hjust = 1, margin = margin(r = 4)),
          axis.text.x      = element_text(family = "Arial", size = 30, face = "plain", colour = "black", angle = 90, hjust = 1, vjust = 1, margin = margin(t = 3)),
          axis.ticks       = element_blank(),
          axis.line        = element_blank(),
          panel.border     = element_rect(fill = NA, colour = "black", linewidth = 2.0),
          panel.grid       = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA),
          plot.title       = element_text(family = "Arial", size = 50, face = "bold.italic", colour = "black", hjust = 0, margin = margin(b = 2)),
          legend.position  = "right",
          legend.title     = element_text(family = "Arial", size = 20, hjust = 0.5),
          legend.text      = element_text(family = "Arial", size = 20),
          legend.margin    = margin(l = 2),
          plot.margin      = margin(t = 16, r = 12, b = 20, l = 12)) +
    ggtitle(full_title)
  ggsave(filename = out_path, plot = p,
         width = 1.6 + n_genes * cell_width + 1.4, height = fig_height,
         dpi = dpi, bg = "white")
  cat(sprintf("  Saved: %s\n", out_path))
  invisible(p)
}

tier1_dir <- file.path(OUT_DIR, "Tier1_Full")
mats_g <- build_lfc_matrix(all_results, glycolysis_genes, tissues)
cap_g  <- compute_lfc_cap(mats_g$lfc, mats_g$padj)
write.csv(as.data.frame(mats_g$lfc) %>% rownames_to_column("cell_line"),
          file.path(tier1_dir, "A_Glycolysis_Full_LFC_Table.csv"), row.names = FALSE)
draw_pathway_heatmap(lfc_mat = mats_g$lfc, padj_mat = mats_g$padj,
                     panel_title = "Glycolysis (Warburg)",
                     subtitle    = paste0("Full gene set  |  n = ", ncol(mats_g$lfc), " genes"),
                     panel_color = "#8B0000",
                     out_path    = file.path(tier1_dir, "Full_A_Glycolysis_Heatmap.png"),
                     lfc_cap     = cap_g)
mats_h <- build_lfc_matrix(all_results, hypoxia_genes, tissues)
cap_h  <- compute_lfc_cap(mats_h$lfc, mats_h$padj)
write.csv(as.data.frame(mats_h$lfc) %>% rownames_to_column("cell_line"),
          file.path(tier1_dir, "B_Hypoxia_Full_LFC_Table.csv"), row.names = FALSE)
draw_pathway_heatmap(lfc_mat = mats_h$lfc, padj_mat = mats_h$padj,
                     panel_title = "Hypoxia Response",
                     subtitle    = paste0("Full gene set  |  n = ", ncol(mats_h$lfc), " genes"),
                     panel_color = "#1F4E79",
                     out_path    = file.path(tier1_dir, "Full_B_Hypoxia_Heatmap.png"),
                     lfc_cap     = cap_h)
mats_c <- build_lfc_matrix(all_results, ctbp_genes, tissues)
cap_c  <- compute_lfc_cap(mats_c$lfc, mats_c$padj)
write.csv(as.data.frame(mats_c$lfc) %>% rownames_to_column("cell_line"),
          file.path(tier1_dir, "C_CtBP_Full_LFC_Table.csv"), row.names = FALSE)
draw_pathway_heatmap(lfc_mat = mats_c$lfc, padj_mat = mats_c$padj,
                     panel_title = "CtBP-related",
                     subtitle    = paste0("Full gene set  |  n = ", ncol(mats_c$lfc), " genes"),
                     panel_color = "#B8860B",
                     out_path    = file.path(tier1_dir, "Full_C_CtBP_Heatmap.png"),
                     lfc_cap     = cap_c)
cat(sprintf("  (padj < %.2f  &  |log2FC| >= %.1f  in >= 1 cell line)\n", PADJ_THRESH, LFC_THRESH))
tier2_dir <- file.path(OUT_DIR, "Tier2_Significant")
sig_g <- filter_sig_genes(mats_g$lfc, mats_g$padj)
if (!is.null(sig_g))
{
  write.csv(as.data.frame(sig_g$lfc) %>% rownames_to_column("cell_line"),
            file.path(tier2_dir, "A_Glycolysis_Sig_LFC_Table.csv"), row.names = FALSE)
  draw_pathway_heatmap(lfc_mat = sig_g$lfc, padj_mat = sig_g$padj,
                       panel_title = "Glycolysis (Warburg)",
                       subtitle    = paste0("Significant DEGs  |  n = ", ncol(sig_g$lfc),
                                            " genes  |  padj < 0.05  &  |log2FC| >= 1"),
                       panel_color = "#8B0000",
                       out_path    = file.path(tier2_dir, "Sig_A_Glycolysis_Heatmap.png"),
                       lfc_cap     = cap_g)
}
sig_h <- filter_sig_genes(mats_h$lfc, mats_h$padj)
if (!is.null(sig_h))
{
  write.csv(as.data.frame(sig_h$lfc) %>% rownames_to_column("cell_line"),
            file.path(tier2_dir, "B_Hypoxia_Sig_LFC_Table.csv"), row.names = FALSE)
  draw_pathway_heatmap(lfc_mat = sig_h$lfc, padj_mat = sig_h$padj,
                       panel_title = "Hypoxia Response",
                       subtitle    = paste0("Significant DEGs  |  n = ", ncol(sig_h$lfc),
                                            " genes  |  padj < 0.05  &  |log2FC| >= 1"),
                       panel_color = "#1F4E79",
                       out_path    = file.path(tier2_dir, "Sig_B_Hypoxia_Heatmap.png"),
                       lfc_cap     = cap_h)
}
sig_c <- filter_sig_genes(mats_c$lfc, mats_c$padj)
if (!is.null(sig_c))
{
  write.csv(as.data.frame(sig_c$lfc) %>% rownames_to_column("cell_line"),
            file.path(tier2_dir, "C_CtBP_Sig_LFC_Table.csv"), row.names = FALSE)
  draw_pathway_heatmap(lfc_mat = sig_c$lfc, padj_mat = sig_c$padj,
                       panel_title = "CtBP-related",
                       subtitle    = paste0("Significant DEGs  |  n = ", ncol(sig_c$lfc),
                                            " genes  |  padj < 0.05  &  |log2FC| >= 1"),
                       panel_color = "#B8860B",
                       out_path    = file.path(tier2_dir, "Sig_C_CtBP_Heatmap.png"),
                       lfc_cap     = cap_c)
}
pathway_summary <- bind_rows(
  all_combined %>% filter(gene %in% glycolysis_genes) %>%
    mutate(pathway = "Aerobic_Glycolysis"),
  all_combined %>% filter(gene %in% hypoxia_genes, !gene %in% glycolysis_genes) %>%
    mutate(pathway = "Hypoxia_Response"),
  all_combined %>% filter(gene %in% ctbp_genes, !gene %in% glycolysis_genes,
                          !gene %in% hypoxia_genes) %>%
    mutate(pathway = "CtBP_related")) %>%
  select(pathway, gene, cell_line, baseMean, log2FoldChange, pvalue, padj, direction) %>%
  arrange(pathway, gene, cell_line)
write.csv(pathway_summary, file.path(OUT_DIR, "Pathway_Genes_DESeq2_Summary.csv"), row.names = FALSE)
cat("  Saved: Pathway_Genes_DESeq2_Summary.csv\n")