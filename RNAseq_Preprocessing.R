# RNA seq preprocessing
# Count matrix, Normalization and Quality control
rm(list = ls())
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(showtext)
font_add("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_auto()
# User Inputs
work_dir <- readline(prompt = "Enter working directory: ")
setwd(work_dir)
result_folder <- "Results"
if(!dir.exists(result_folder)) dir.create(result_folder)
min_total_count <- 10
# Data inputs
tissues <- c("A549", "HCT116", "PC3", "T47D")
conditions <- c("Nx", "Hx")
n_rep <- 3
# Locate and check the count files
esc <- function(x) gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", x)
all_txt <- list.files(".", pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
find_files <- function(tissue, cond)
{
  pat <- paste0("^", esc(tissue), "_", esc(cond), "([0-9]+)\\.txt$")
  hits <- all_txt[grepl(pat, basename(all_txt))]
  if (length(hits) == 0) return(NULL)
  data.frame(tissue = tissue, condition = cond, replicate = as.integer(sub(pat, "\\1", basename(hits))), file = hits, stringsAsFactors = FALSE)
}
found <- list()
for (t in tissues) for (cn in conditions) found[[length(found) + 1]] <- find_files(t, cn)
meta <- do.call(rbind, found)
if (is.null(meta) || nrow(meta) == 0)
{
  stop("No files named <tissue>_<condition><replicate>.txt found", getwd(), "\nExpected format")
}
missing_combo <- setdiff(paste(rep(tissues, each = length(conditions)), conditions, sep = "/"), paste(meta$tissue, meta$condition, sep = "/"))
if (length(missing_combo) > 0)
{
  stop("Fill missing for: ", paste(missing_combo, collapse = ", "))
}
for (t in tissues)
{
  for (cn in conditions)
  {
    reps_found <- sort(meta$replicate[meta$tissue == t & meta$condition == cn])
    expected <- seq_len(n_rep)
    if (!identical(reps_found, expected))
    {
      stop(sprintf("%s_%s: expected replicates %s, found replicates %s", t, cn, paste(expected, collapse = ","), paste(reps_found, collapse = ",")))
    } }}
if (anyDuplicated(meta[, c("tissue", "condition", "replicate")]) > 0)
{
  stop("Duplicate files found")
}

meta$sample <- paste0(meta$tissue, "_", meta$condition, meta$replicate)
meta$sampleintissue <- paste0(meta$condition, meta$replicate)
meta <- meta[order(match(meta$tissue, tissues), match(meta$condition, conditions), meta$replicate), ]
rownames(meta) <- NULL
cat("\nReplicates found per tissue x condition:\n")
print(table(tissue    = factor(meta$tissue,     levels = tissues), condition = factor(meta$condition,  levels = conditions)))

# Read the count files
read_counts <- function(file)
{
  x <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 2 || !is.numeric(x[[2]]))
    stop("Unexpected format in ", file)
  x <- x[!grepl("^__", x[[1]]), 1:2]
  if (anyDuplicated(x[[1]])) stop("Duplicated gene IDs in ", file)
  x
}
raw <- setNames(lapply(meta$file, read_counts), meta$sample)
# Merge counts
merge_counts <- function(ids, new_names)
{
  lst <- Map(function(id, nm)
  {
    x <- raw[[id]]; colnames(x) <- c("Gene", nm); x
  }, ids, new_names)
  m <- Reduce(function(x, y) merge(x, y, by = "Gene"), lst)
  if (nrow(m) < max(vapply(raw[ids], nrow, integer(1))))
    warning("Genes absent from some files were dropped: ", nrow(m), " genes are shared by all samples", call. = FALSE)
  m
}

# Processing data
process_dataset <- function(count_matrix, col_data, design, out_dir, prefix, title, pca_color, pca_shape = NULL, pca_label_by = NULL)
{
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Count matrix
  cm_file <- file.path(out_dir, paste0(prefix, "_Count_Matrix.csv"))
  write.csv(count_matrix, cm_file, row.names = FALSE)
  cat("Saved:", cm_file, "\n")
  # Normalization
  cnt <- as.matrix(count_matrix[, -1])
  rownames(cnt) <- count_matrix$Gene
  col_data <- col_data[colnames(cnt), , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(countData = cnt, colData = col_data, design = design)
  dds <- dds[rowSums(counts(dds)) >= min_total_count, ]
  dds <- estimateSizeFactors(dds)
  norm_file <- file.path(out_dir, paste0(prefix, "_Normalized_Counts.csv"))
  write.csv(counts(dds, normalized = TRUE), norm_file)
  saveRDS(dds, file.path(out_dir, paste0(prefix, "_dds.rds")))
  cat("Saved:", norm_file, "\n")
  # Variance-stabilizing transformation
  vsd <- if (nrow(dds) >= 1000) 
  {
    vst(dds, blind = TRUE)
  } 
  else 
  {
    varianceStabilizingTransformation(dds, blind = TRUE)
  }
  # QC1 - Library size
  png(file.path(out_dir, paste0(prefix, "_QC1_Library_Size.png")), width = 2200, height = 1800, res = 300)
  par(mar = c(5, 12, 5, 2), family = "Arial")
  barplot(colSums(cnt), horiz = TRUE, las = 1, col = "steelblue", xlab = "Raw Counts", main = paste(title, "Library Size"))
  dev.off()
  # QC2 - PCA
  intgroup_vars <- unique(c(pca_color, pca_shape, pca_label_by))
  intgroup_vars <- intgroup_vars[!sapply(intgroup_vars, is.null)]
  pca <- plotPCA(vsd, intgroup = intgroup_vars, returnData = TRUE)
  pv  <- round(100 * attr(pca, "percentVar"))
  p <- ggplot(pca, aes(PC1, PC2, colour = .data[[pca_color]])) + geom_point(size = 3) + xlab(paste0("PC1: ", pv[1], "% variance")) + ylab(paste0("PC2: ", pv[2], "% variance")) + ggtitle(paste(title, "PCA")) + theme_bw() +
    theme(text = element_text(family = "Arial",),
          axis.title = element_text(family = "Arial", size = 32),
          axis.text = element_text(family = "Arial", size = 30),
          legend.text = element_text(family = "Arial", size = 30),
          legend.title = element_text(family = "Arial", size = 30),
          plot.title = element_text(family = "Arial", size = 35))
  if (!is.null(pca_shape)) p <- p + aes(shape = .data[[pca_shape]])
  if (!is.null(pca_label_by)) 
  {
    # Combined PCA 
    centroids <- aggregate(cbind(PC1, PC2) ~ ., data = pca[, c(pca_label_by, "PC1", "PC2")], FUN = mean)
    colnames(centroids)[1] <- pca_label_by
    p <- p + geom_text(data = centroids, aes(x = PC1, y = PC2, label = .data[[pca_label_by]]),
                       inherit.aes = FALSE, colour = "black", fontface = "bold", size = 10, vjust = -1.2, family = "Arial")
  } 
  else 
  {
    # Tissue-wise PCA
    p <- p + geom_text(aes(label = name), vjust = -1, size = 10, show.legend = FALSE, family = "Arial")
  }
  ggsave(file.path(out_dir, paste0(prefix, "_QC2_PCA.png")), p, width = 6.5, height = 5, dpi = 300)
  # QC3 - Sample distance heatmap
  sampleDist <- dist(t(assay(vsd)))
  png(file.path(out_dir, paste0(prefix, "_QC3_Sample_Distance_Heatmap.png")), width = 1800, height = 1600, res = 300)
  pheatmap(as.matrix(sampleDist), main = paste(title, "Sample Distance"), fontfamily = "Arial")
  dev.off()
  # QC4 - Correlation heatmap
  corMatrix <- cor(assay(vsd))
  png(file.path(out_dir, paste0(prefix, "_QC4_Correlation_Heatmap.png")), width = 1800, height = 1600, res = 300)
  pheatmap(corMatrix, main = paste(title, "Sample Correlation"), fontfamily = "Arial")
  dev.off()
  
  invisible(dds)
}

dir.create(result_folder, showWarnings = FALSE)

# Tissue-wise analysis
for (tissue in tissues)
{
  cat("\nTissue-wise analysis:", tissue, "\n")
  m_t <- meta[meta$tissue == tissue, ]
  cm  <- merge_counts(m_t$sample, m_t$sampleintissue)
  cd  <- data.frame(condition = factor(m_t$condition, levels = conditions), replicate = factor(m_t$replicate), row.names = m_t$sampleintissue)
  process_dataset(cm, cd, ~ condition, file.path(result_folder, tissue), tissue, tissue, pca_color = "condition")
}
# Combined analysis
cat("\nCombined analysis\n")
cm_all <- merge_counts(meta$sample, meta$sample)
cd_all <- data.frame(tissue = factor(meta$tissue, levels = tissues), condition = factor(meta$condition, levels = conditions), replicate = factor(meta$replicate), row.names = meta$sample)
process_dataset(cm_all, cd_all, ~ tissue + condition, file.path(result_folder, "Combined"), "Combined", "Combined", pca_color = "condition", pca_label_by = "tissue")      

# Save settings for the DESeq2 analysis
saveRDS(list(work_dir = getwd(),
             result_folder = result_folder,
             tissues = tissues,
             conditions = conditions,
             min_total_count = min_total_count,
             metadata  = meta),
        file.path(result_folder, "analysis_config.rds")
)
write.csv(meta, file.path(result_folder, "Sample_Metadata.csv"), row.names = FALSE)
cat("\nPreprocessing Completed")