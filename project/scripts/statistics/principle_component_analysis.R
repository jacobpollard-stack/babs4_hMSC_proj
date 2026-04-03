# ==============================================================================
# Livecyte Data Analysis
# - Statistical analysis:
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
# 3a.Select numeric columns for PCA and remove replicate and tracking.id
#
numeric_metrics <- metrics |> 
  select(-replicate, -tracking.id) |> 
  select(where(is.numeric))
#
# 3b. Perform PCA
#
pca_result <- prcomp(numeric_metrics, scale. = TRUE)
#
# Create a data frame for PCA results
#
pca_df <- as.data.frame(pca_result$x) |> 
  bind_cols(metrics |> select(-where(is.numeric)))
#
# 3c. Plot PCA results in 2d with ellipses
#
ggplot(pca_df, aes(x = PC1, y = PC2, colour = clone)) +
  geom_point(size = 3) +
  labs(x = 'Morphology',
       y = 'Movement') +
  theme_minimal() +
  stat_ellipse(level = 0.95, aes(fill = clone), alpha = 0, geom = "polygon", show.legend = FALSE) +
  theme(legend.position = "right")
#
# We see that dry.mass and volume are perfectly correlated. This is expected as dry.mass is a function of volume, calculated by the microscope.
#
# Plot 3d plot of PC1, PC2, and PC3, representing morphology, movement, and persistence, respectively.
#
plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~clone, type = 'scatter3d', mode = 'markers', marker = list(size = 5, opacity = 0.75)) |> 
  layout(scene = list(xaxis = list(title = 'Morphology'),
                      yaxis = list(title = 'Movement'),
                      zaxis = list(title = 'Persistence')))
#
# The 3D plot shows that the clones form two seperate clouds in the PCA space, with some overlap. This indicates that, when considering morphology, movement, and 'persistence' (how long it stays on screen, and its mean speed), the two clones are largely distinct. This suggests that the clones have different phenotypic profiles, which could be due to different osteogenic potential.