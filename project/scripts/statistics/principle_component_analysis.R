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
# 3c. Plot PCA results in 2d
#
ggplot(pca_df, aes(x = PC1, y = PC2, colour = clone)) +
  geom_point(size = 3) +
  labs(x = 'Morphology',
       y = 'Movement')
#
# We see that dry.mass and volume are perfectly correlated. This is expected as dry.mass is a function of volume, calculated by the microscope.
#
# Plot 3d plot of PC1, PC2, and PC3, representing morphology, movement, and persistence, respectively.
#
plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~clone, type = 'scatter3d', mode = 'markers', marker = list(size = 5, opacity = 0.75)) |> 
  layout(scene = list(xaxis = list(title = 'Morphology'),
                      yaxis = list(title = 'Movement'),
                      zaxis = list(title = 'Persistence')))
