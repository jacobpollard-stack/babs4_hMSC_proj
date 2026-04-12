# ==============================================================================
# Livecyte Data Analysis
# - Statistical analysis and figures:
# - Principle Component Analysis (PCA)
# ==============================================================================
#
# 1. Load libraries
#
library(tidyverse)
library(ggplot2)
library(plotly)
#
# 2. Load data
#
metrics <- read_tsv('project/data/movement_morphology/livecyte_collapsed_filtered.tsv')
#
# 3. Principle Component Analysis (PCA)
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
# 3d. Perform statistical tests to determine which PCs explain the most variance in the data
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
# We see that PC1, 2, and 3 explain ~78 of the variance in the data, so we will focus on these three PCs for further analysis. PCs after this become less biologically interpretable
#
# 4. Peform non-parametric statistical tests to determine the difference in Principle Components between the two clones
#
# 4a. Isolate the 3 most biologically interpretable PCs: PC1 largely represents morphology; PC2 mostly represents how far the cell moved and how much it meandered, but also seems to be correlated with length to width ratio and dry mass. As dry mass, volume, and total path length are not significantly different between the two clones, PC2 could be identifying debris: A high PC2 value could be a debris marker. PC3 largely represents movement, which is also correlated with other morphological parameters.
#
PC1_PC2_PC3 <- pca_df |> 
  select(clone, PC1, PC2, PC3)
#
# 4b. Perform Wilcoxon rank-sum test for each PC
#
wilcox_PC1 <- wilcox.test(PC1 ~ clone, data = PC1_PC2_PC3)
wilcox_PC2 <- wilcox.test(PC2 ~ clone, data = PC1_PC2_PC3)
wilcox_PC3 <- wilcox.test(PC3 ~ clone, data = PC1_PC2_PC3)
#
# 4c. Extract median values for each clone for each PC
med_PC2_A <- median(pca_df$PC2[pca_df$clone == "A"])
med_PC2_B <- median(pca_df$PC2[pca_df$clone == "B"])
med_PC1_A <- median(pca_df$PC1[pca_df$clone == "A"])
med_PC1_B <- median(pca_df$PC1[pca_df$clone == "B"])
med_PC3_A <- median(pca_df$PC3[pca_df$clone == "A"])
med_PC3_B <- median(pca_df$PC3[pca_df$clone == "B"])
#
# 5. Plot results with p-values
#
# 5a. PC2 on PC1
#
# 5ai. Build base plot with median points for each clone
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
# 5aii. Statistical labels
#
  annotate("segment", x = -4.5, xend = -4.5, y = med_PC2_A, yend = med_PC2_B) +
  annotate("segment", x = -4.5, xend = -4.5 + 0.1, y = med_PC2_A, yend = med_PC2_A) +
  annotate("segment", x = -4.5, xend = -4.5 + 0.1, y = med_PC2_B, yend = med_PC2_B) +
  annotate("text", x = -4.8, y = -0.1, label = "***", size = 4.5) +
#
  annotate("segment", x = med_PC1_B, xend = med_PC1_A, y = -3.5) +
  annotate("segment", x = med_PC1_B, y = -3.5, yend = -3.5 + 0.35) +
  annotate("segment", x = med_PC1_A, y = -3.5, yend = -3.5 + 0.35) +
  annotate("text", x = 0.1, y = -4, label = "***", size = 4.5)
#
# 5b. PC3 on PC1
#
# 5bi. Build base plot with median points for each clone
#
PC3_PC1_plot <- ggplot(pca_df, aes(x = PC1, y = PC3, colour = clone, alpha = clone)) +
  geom_point(size = 3) +
  scale_alpha_manual(values = c("A" = 0.4, "B" = 0.25), guide = "none") +
  labs(x = 'PC1',
       y = 'PC3',
       colour = "Clone") +
  theme_test() +
  stat_ellipse(level = 0.95, aes(fill = clone), alpha = 0, geom = "polygon", show.legend = FALSE) +
  theme(legend.position = 'none') +
  # Median point for clone A
  geom_point(aes(x = med_PC1_A, y = med_PC3_A), fill = "red", colour = "red", shape = 24, size = 3) +
  # Median point for clone B
  geom_point(aes(x = med_PC1_B, y = med_PC3_B), colour = "blue", fill = "blue", shape = 24, size = 3) +
#
# 5bii. Statistical labels
#
  annotate("segment", x = -4.5, xend = -4.5, y = med_PC3_A, yend = med_PC3_B) +
  annotate("segment", x = -4.5, xend = -4.5 + 0.1, y = med_PC3_A, yend = med_PC3_A) +
  annotate("segment", x = -4.5, xend = -4.5 + 0.1, y = med_PC3_B, yend = med_PC3_B) +
  annotate("text", x = -4.8, y = 0, label = "***", size = 4.5) +
  #
  annotate("segment", x = med_PC1_B, xend = med_PC1_A, y = -3) +
  annotate("segment", x = med_PC1_B, y = -3, yend = -3 + 0.2) +
  annotate("segment", x = med_PC1_A, y = -3, yend = -3 + 0.2) +
  annotate("text", x = 0.1, y = -3.3, label = "***", size = 4.5)
#
# 5c. PC3 on PC2
#
# 5ci. Build base plot with median points for each clone
#
PC3_PC2_plot <- ggplot(pca_df, aes(x = PC2, y = PC3, colour = clone, alpha = clone)) +
  geom_point(size = 3) +
  scale_alpha_manual(values = c("A" = 0.4, "B" = 0.25), guide = "none") +
  labs(x = 'PC2',
       y = 'PC3',
       colour = "Clone") +
  theme_test() +
  stat_ellipse(level = 0.95, aes(fill = clone), alpha = 0, geom = "polygon", show.legend = FALSE) +
  theme(legend.position = 'none') +
  # Median point for clone A
  geom_point(aes(x = med_PC2_A, y = med_PC3_A), fill = "red", colour = "red", shape = 24, size = 3) +
  # Median point for clone B
  geom_point(aes(x = med_PC2_B, y = med_PC3_B), colour = "blue", fill = "blue", shape = 24, size = 3) +
#
# 5bii. Statistical labels
#
  annotate("segment", x = -3.2, xend = -3.2, y = med_PC3_A, yend = med_PC3_B) +
  annotate("segment", x = -3.2, xend = -3.2 + 0.1, y = med_PC3_A, yend = med_PC3_A) +
  annotate("segment", x = -3.2, xend = -3.2 + 0.1, y = med_PC3_B, yend = med_PC3_B) +
  annotate("text", x = -3.5, y = 0, label = "***", size = 4.5) +
  #
  annotate("segment", x = med_PC2_B, xend = med_PC2_A, y = -3) +
  annotate("segment", x = med_PC2_B, y = -3, yend = -3 + 0.2) +
  annotate("segment", x = med_PC2_A, y = -3, yend = -3 + 0.2) +
  annotate("text", x = -0.15, y = -3.3, label = "***", size = 4.5)
#
# 6. Build a loadings plot to show which features contribute most to each PC
#
# 6a. Extract the loadings for the first 3 PCs
#
loadings_mat <- pca_result$rotation[, 1:3]
#
# 6b. Convert loadings to long format
#
loadings_long <- as.data.frame(loadings_mat) |>
  rownames_to_column("feature") |>
  pivot_longer(cols = PC1:PC3, names_to = "PC", values_to = "loading") |>
  mutate(feature = factor(feature, levels = rev(rownames(loadings_mat))))
#
# 6c. Plot loadings as a heatmap
#
loadings_heatmap <- ggplot(loadings_long, aes(x = PC, y = feature, fill = loading)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(loading, 2)), size = 3.2) +
  scale_fill_gradient2(low = "#ED7117", mid = "white", high = "#6F2DA8",
                       midpoint = 0, limits = c(-0.6, 0.6)) +
  labs(x = NULL, y = NULL, fill = "Loading") +
  theme_test() +
  theme(legend.position = "right")
#
# 7. Export all plots
#
ggsave('project/figures/PCA/pca_PC2_PC1_plot.png', PC2_PC1_plot, width = 5, height = 4)
ggsave('project/figures/PCA/pca_PC3_PC1_plot.png', PC3_PC1_plot, width = 5, height = 4)
ggsave('project/figures/PCA/pca_PC3_PC2_plot.png', PC3_PC2_plot, width = 5, height = 4)
ggsave('project/figures/PCA/pca_loadings_heatmap.png', loadings_heatmap, width = 5, height = 3.75)