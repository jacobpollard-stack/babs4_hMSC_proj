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
# 3a. Select numeric columns for PCA and remove replicate and tracking.id
#
numeric_metrics <- metrics |> 
  select(-replicate, -tracking.id) |> 
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
# We see that PC1, 2, and 3 explain ~70% of the variance in the data, so we will focus on these three PCs for further analysis. PCs after this become less interpretable
#
# 4. Peform non-parametric statistical tests to determine the difference in Principle Components between the two clones
#
# 4a. Isolate the 3 most biologically interpretable PCs: PC1 largely represents morphology; PC2 largely represents movement; PC3 largely represents 'persistence' (how long the cell stays on screen, and its mean speed).
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
  labs(x = 'PC1: Morphology',
       y = 'PC2: Movement',
       colour = "Clone") +
  theme_test() +
  stat_ellipse(level = 0.95, aes(fill = clone), alpha = 0, geom = "polygon", show.legend = FALSE) +
  theme(legend.position = "right") +
# Median point for clone A
  geom_point(aes(x = med_PC1_A, y = med_PC2_A), fill = "red", colour = "red", shape = 24, size = 3) +
# Median point for clone B
  geom_point(aes(x = med_PC1_B, y = med_PC2_B), colour = "blue", fill = "blue", shape = 24, size = 3)
#
# 5aii. Statistical labels
#
PC2_PC1_plot + 
  annotate("segment", x = -4.5, xend = -4.5, y = med_PC2_A, yend = med_PC2_B) +
  annotate("segment", x = -4.5, xend = -4.5 + 0.1, y = med_PC2_A, yend = med_PC2_A) +
  annotate("segment", x = -4.5, xend = -4.5 + 0.1, y = med_PC2_B, yend = med_PC2_B) +
  annotate("text", x = -4.8, y = 0.35, label = "***", size = 4.5) +
#
  annotate("segment", x = med_PC1_B, xend = med_PC1_A, y = -4.2) +
  annotate("segment", x = med_PC1_B, y = -4.2, yend = -4.2 + 0.35) +
  annotate("segment", x = med_PC1_A, y = -4.2, yend = -4.2 + 0.35) +
  annotate("text", x = 0, y = -4.6, label = "***", size = 4.5)
#
# 5b. PC3 on PC1
#
# 5bi. Build base plot with median points for each clone
#
PC3_PC1_plot <- ggplot(pca_df, aes(x = PC1, y = PC3, colour = clone, alpha = clone)) +
  geom_point(size = 3) +
  scale_alpha_manual(values = c("A" = 0.4, "B" = 0.25), guide = "none") +
  labs(x = 'PC1: Morphology',
       y = 'PC3: Persistence',
       colour = "Clone") +
  theme_test() +
  stat_ellipse(level = 0.95, aes(fill = clone), alpha = 0, geom = "polygon", show.legend = FALSE) +
  theme(legend.position = "right") +
  # Median point for clone A
  geom_point(aes(x = med_PC1_A, y = med_PC3_A), fill = "red", colour = "red", shape = 24, size = 3) +
  # Median point for clone B
  geom_point(aes(x = med_PC1_B, y = med_PC3_B), colour = "blue", fill = "blue", shape = 24, size = 3)
#
# 5bii. Statistical labels
#
PC3_PC1_plot + 
  annotate("segment", x = -4.5, xend = -4.5, y = med_PC3_A, yend = med_PC3_B) +
  annotate("segment", x = -4.5, xend = -4.5 + 0.1, y = med_PC3_A, yend = med_PC3_A) +
  annotate("segment", x = -4.5, xend = -4.5 + 0.1, y = med_PC3_B, yend = med_PC3_B) +
  annotate("text", x = -4.75, y = -0.2, label = "***", size = 4.5) +
  #
  annotate("segment", x = med_PC1_B, xend = med_PC1_A, y = -3.5) +
  annotate("segment", x = med_PC1_B, y = -3.5, yend = -3.5 + 0.2) +
  annotate("segment", x = med_PC1_A, y = -3.5, yend = -3.5 + 0.2) +
  annotate("text", x = 0, y = -3.8, label = "***", size = 4.5)
#
# 5c. PC3 on PC2
#
# 5ci. Build base plot with median points for each clone
#
PC3_PC2_plot <- ggplot(pca_df, aes(x = PC2, y = PC3, colour = clone, alpha = clone)) +
  geom_point(size = 3) +
  scale_alpha_manual(values = c("A" = 0.4, "B" = 0.25), guide = "none") +
  labs(x = 'PC2: Movement',
       y = 'PC3: Persistence',
       colour = "Clone") +
  theme_test() +
  stat_ellipse(level = 0.95, aes(fill = clone), alpha = 0, geom = "polygon", show.legend = FALSE) +
  theme(legend.position = "right") +
  # Median point for clone A
  geom_point(aes(x = med_PC2_A, y = med_PC3_A), fill = "red", colour = "red", shape = 24, size = 3) +
  # Median point for clone B
  geom_point(aes(x = med_PC2_B, y = med_PC3_B), colour = "blue", fill = "blue", shape = 24, size = 3)
#
# 5bii. Statistical labels
#
PC3_PC2_plot + 
  annotate("segment", x = -4.3, xend = -4.3, y = med_PC3_A, yend = med_PC3_B) +
  annotate("segment", x = -4.3, xend = -4.3 + 0.1, y = med_PC3_A, yend = med_PC3_A) +
  annotate("segment", x = -4.3, xend = -4.3 + 0.1, y = med_PC3_B, yend = med_PC3_B) +
  annotate("text", x = -4.6, y = -0.2, label = "***", size = 4.5) +
  #
  annotate("segment", x = med_PC2_B, xend = med_PC2_A, y = -3.5) +
  annotate("segment", x = med_PC2_B, y = -3.5, yend = -3.5 + 0.2) +
  annotate("segment", x = med_PC2_A, y = -3.5, yend = -3.5 + 0.2) +
  annotate("text", x = 0.375, y = -3.8, label = "***", size = 4.5)
#
# 5. Build a loadings plot to show which features contribute most to each PC
#
# 5a. Extract the loadings for the first 3 PCs
#
loadings_mat <- pca_result$rotation[, 1:3]
#
# 5b. Create a heatmap of the loadings
#
rownames(loadings_mat) <- c("dry.mass", "volume", "radius", "sphericity", "length.to.width.ratio", "mean.speed", "total_path_length", "final_displacement")
#
# 5c. Convert loadings to long format
#
loadings_long <- as.data.frame(loadings_mat) |>
  rownames_to_column("feature") |>
  pivot_longer(cols = PC1:PC3, names_to = "PC", values_to = "loading") |>
  mutate(feature = factor(feature, levels = rev(rownames(loadings_mat))))

loadings_heatmap <- ggplot(loadings_long, aes(x = PC, y = feature, fill = loading)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(loading, 2)), size = 3.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-0.6, 0.6)) +
  labs(x = NULL, y = NULL, fill = "Loading") +
  theme_test() +
  theme(legend.position = "right")











# 3d. Plot results
#
# 3di. PC2 on PC1
#
ggplot(pca_df, aes(x = PC1, y = PC2, colour = clone)) +
  geom_point(size = 3) +
  labs(x = 'PC1: Morphology',
       y = 'PC2: Movement',
       colour = "Clone") +
  theme_test() +
  stat_ellipse(level = 0.95, aes(fill = clone), alpha = 0, geom = "polygon", show.legend = FALSE) +
  theme(legend.position = "right")
#
# 3dii. PC3 on PC1
#
ggplot(pca_df, aes(x = PC1, y = PC3, colour = clone)) +
  geom_point(size = 3) +
  labs(x = 'PC1: Morphology',
       y = 'PC3: Persistence',
       colour = "Clone") +
  theme_test() +
  stat_ellipse(level = 0.95, aes(fill = clone), alpha = 0, geom = "polygon", show.legend = FALSE) +
  theme(legend.position = "right")
#
# 3diii. PC3 on PC2
#
ggplot(pca_df, aes(x = PC2, y = PC3, colour = clone)) +
  geom_point(size = 3) +
  labs(x = 'PC2: Movement',
       y = 'PC3: Persistence',
       colour = "Clone") +
  theme_test() +
  stat_ellipse(level = 0.95, aes(fill = clone), alpha = 0, geom = "polygon", show.legend = FALSE) +
  theme(legend.position = "right")
#


# Plot 3d plot of PC1, PC2, and PC3, representing morphology, movement, and persistence, respectively
#
plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~clone, colors = c('blue', 'red'), type = 'scatter3d', mode = 'markers', marker = list(size = 5), opacity = 0.5) |> 
  layout(scene = list(xaxis = list(title = 'PC1: Morphology'),
                      yaxis = list(title = 'PC2: Movement'),
                      zaxis = list(title = 'PC3: Persistence')),
         legend = list(title = list(text = 'Clone')))
#
# The 3D plot shows that the clones form two seperate clouds in the PCA space, with some overlap. This indicates that, when considering morphology, movement, and 'persistence' (how long it stays on screen, and its mean speed), the two clones are largely distinct. This suggests that the clones have different phenotypic profiles, which could be due to different osteogenic potential.
#
