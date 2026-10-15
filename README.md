# Paralog-specific bistability links oxygen tension to CtBP-driven transcriptional reprogramming in tumor cells
## Saumya Shukla, Raghvendra Singh

## Requirements
•	R (v4.0 or higher recommended)  
•	Libraries: DESeq2, ggplot2, pheatmap, dplyr, tidyr, tibble, scales, ashr, showtext  
•	MATLAB (R2025a recommended)  
•	No additional toolboxes needed, uses built-in ode15s

## Files
### aerobic_glycolysis.m
Kinetic ODE model of aerobic glycolysis simulating 21 metabolic species across six oxygen levels, from near-anoxia to hyperoxia. Tracks how cytosolic NADH changes at steady state. Prints steady-state NADH values for each oxygen level and bar chart of total and free cytosolic NADH across oxygen conditions. 
### CtBP_oligomerization.m
Models CtBP1 and CtBP2 tetramer assembly as a function of free cytosolic NADH using Hill equations. Uses NADH values from the aerobic glycolysis model across six oxygen states. Prints tetramer fraction (%) for CtBP1 and CtBP2 at each oxygen state, grouped bar chart comparing tetramer fractions (with uncertainty range bars) and Hill curve plots for CtBP1 and CtBP2 showing tetramer fraction vs. free NADH.
### CtBP_utilization.m
Compares three models of CtBP utilization as a function of free NADH: linear, Michaelis-Menten, and cooperative (Hill) against the Hill-equation production curve. Done separately for CtBP1 and CtBP2. Gives figures (one per paralog) showing all three utilization models alongside the actual production curve, with oxygen state data points annotated. 
### RNAseq_Preprocessing.R
Reads raw RNA-seq count files, builds a count matrix, normalizes the data, and runs quality control. Gives count matrix and normalized counts (.csv), QC plots: library size bar chart, PCA plot, sample distance heatmap, correlation heatmap (.png) and DESeq2 object and analysis config (.rds) for use in the next script.  
Data needed: RNA seq data GSE200204 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200024)
### RNAseq_analysis.R
Runs DESeq2 differential gene expression analysis (hypoxia vs. normoxia) across multiple cell lines. Generates pathway-level heatmaps for three gene sets: aerobic glycolysis, hypoxia response, and CtBP-related genes. Produces DESeq2 results per cell line (.csv), full and significant-only heatmaps for each pathway (.png) and summary table of pathway gene results (.csv).  
Requires analysis_config.rds produced by RNAseq_Preprocessing.R
