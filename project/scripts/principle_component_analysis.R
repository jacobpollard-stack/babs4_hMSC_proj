# ====================================================================
# Livecyte Data Analysis
# - Statistical analysis and figures:
# - Principle Component Analysis (PCA)
# ====================================================================
#
# 1. Load libraries --------------------------------------------------
#
library(tidyverse) # for data manipulation and visualization
library(ggplot2) # for plotting
library(plotly) # for interactive plots
library(readxl) # for reading Excel files
#
# 2. Load data -------------------------------------------------------
#
# 2a. Load the filtered livecyte data
#
metrics <- read_tsv('project/data/movement_morphology/livecyte_collapsed_filtered.tsv')
#
# 2b. Load the osteogenesis data and pivot
#
osteo <- read_xlsx("project/data/differentiation/osteogenesis_processed.xlsx") |> 
  pivot_longer(cols = c("0", "8", "8osteo"), names_to = "day", values_to = "absorption") 
#
# 3. Principle Component Analysis (PCA) ------------------------------
#
# 3a. Select numeric columns for PCA and remove replicate and tracking.id, along with parameters that are non-biological
#
numeric_metrics <- metrics |> 
  select(-replicate, -tracking.id, -n_frames, -start_frame) |> 
  select(where(is.numeric))
#
# 3b. Perform PCA
#
pca_result <- prcomp(numeric_metrics, scale. = TRUE)
#
# 3c. Create a data frame for PCA results
#
pca_df <- as.data.frame(pca_result$x) |> 
  bind_cols(metrics |> select(-where(is.numeric)))
#
# 3d. Perform statistical tests to determine which Principle Components (PCs) explain the most variance in the data
#
# 3di. Calculate the proportion of variance explained by each PC
#
pve <- pca_result$sdev^2 / sum(pca_result$sdev^2)
#
# 3dii. Create a scree plot to visualize the proportion of variance explained by each PC
#
data.frame(PC = paste0("PC", 1:length(pve)), PVE = pve) |> 
  ggplot(aes(x = PC, y = PVE)) +
  geom_bar(stat = "identity")
#
# We see that PC1 and 2explain ~60% of the variance in the data, so we will focus on these three PCs for further analysis. PCs after this become less biologically interpretable
#
# 4. Peform non-parametric statistical tests to determine the difference in Principle Components between the two clones ------------
#
# 4a. Identify the 2 most biologically interpretable PCs: PC1 is largely loaded with all the morphological parameters, PC2 seems to also largely be morphology-loaded, but also considers total path length and final displacement. Neither PCs are loaded on mean speed, but this was not significantly different between the clones, so this is not a concern.
#
# 4b. Perform Wilcoxon rank-sum test for each PC
#
wilcox_PC1 <- wilcox.test(PC1 ~ clone, data = pca_df)
wilcox_PC2 <- wilcox.test(PC2 ~ clone, data = pca_df)
#
# 4c. Extract median values for each clone for each PC
#
med_PC2_A <- median(pca_df$PC2[pca_df$clone == "A"])
med_PC2_B <- median(pca_df$PC2[pca_df$clone == "B"])
med_PC1_A <- median(pca_df$PC1[pca_df$clone == "A"])
med_PC1_B <- median(pca_df$PC1[pca_df$clone == "B"])
#
# 5. Plot results with p-values --------------------------------------
#
# 5a. PC2 on PC1: Build base plot with median points for each clone
#
PC2_PC1_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = clone, alpha = clone)) +
  geom_point(size = 3) +
  scale_alpha_manual(values = c("A" = 0.4, "B" = 0.25), guide = "none") +
  labs(x = 'PC1',
       y = 'PC2',
       colour = "Clone") +
  theme_test() +
  stat_ellipse(level = 0.95, aes(fill = clone), alpha = 0, geom = "polygon", show.legend = FALSE) +
  theme(legend.position = 'none') +
# Median point for clone A
  geom_point(aes(x = med_PC1_A, y = med_PC2_A), fill = "red", colour = "red", shape = 24, size = 3) +
# Median point for clone B
  geom_point(aes(x = med_PC1_B, y = med_PC2_B), colour = "blue", fill = "blue", shape = 24, size = 3) +
#
# 5b. Statistical labels
#
  annotate("segment", x = -4.1, xend = -4.1, y = med_PC2_A, yend = med_PC2_B) +
  annotate("segment", x = -4.1, xend = -4.1 + 0.1, y = med_PC2_A, yend = med_PC2_A) +
  annotate("segment", x = -4.1, xend = -4.1 + 0.1, y = med_PC2_B, yend = med_PC2_B) +
  annotate("text", x = -4.4, y = -0.2, label = "***", size = 4.5) +
#
  annotate("segment", x = med_PC1_B, xend = med_PC1_A, y = -3.75) +
  annotate("segment", x = med_PC1_B, y = -3.75, yend = -3.75 + 0.35) +
  annotate("segment", x = med_PC1_A, y = -3.75, yend = -3.75 + 0.35) +
  annotate("text", x = 0, y = -4.2, label = "***", size = 4.5)
PC2_PC1_plot
#
# 6. Build a loadings plot to show which features contribute most to each PC --------------------------------------------------------------
#
# 6a. Extract the loadings for the first 2 PCs
#
loadings_mat <- pca_result$rotation[, 1:2]
#
# 6b. Convert loadings to long format
#
loadings_long <- as.data.frame(loadings_mat) |>
  rownames_to_column("feature") |>
  pivot_longer(cols = PC1:PC2, names_to = "PC", values_to = "loading") |>
  mutate(feature = factor(feature, levels = rev(rownames(loadings_mat))))
#
# 6c. Plot loadings as a heatmap
#
loadings_heatmap <- ggplot(loadings_long, aes(x = PC, y = feature, fill = loading)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(loading, 2)), size = 3.2) +
  scale_fill_gradient2(low = "#ED7117", mid = "white", high = "#6F2DA8",
                       midpoint = 0, limits = c(-0.6, 0.7)) +
  labs(x = NULL, y = NULL, fill = "Loading") +
  theme_test() +
  theme(legend.position = "right")
loadings_heatmap
#
# 7. Export all plots ------------------------------------------------
#
ggsave('project/figures/PCA/pca_PC2_PC1_plot.png', PC2_PC1_plot, width = 5, height = 4)
ggsave('project/figures/PCA/pca_loadings_heatmap.png', loadings_heatmap, width = 5, height = 3.75)
