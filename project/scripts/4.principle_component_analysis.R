# Principal component analysis of livecyte features --------------

# ================================================================

# Experimental overview ------------------------------------------

# We have filtered livecyte data for two clonal cell lines (A and B)
# with nine morphological and movement features. PCA is used to
# identify the main sources of variance in the data and to determine
# which features contribute most to differences between the clones.


# Description of data --------------------------------------------

# The filtered livecyte dataset (livecyte_collapsed_filtered.tsv)
# contains one row per tracking.id with columns: clone, replicate,
# tracking.id, n_frames, start_frame, and nine numeric features.
# The osteogenesis dataset (osteogenesis_processed.xlsx) contains ALP
# absorbance data for comparison. Data is stored in
# project/data/movement_morphology/ and project/data/differentiation/.


# Analysis overview ----------------------------------------------

# This script performs PCA on the filtered livecyte dataset to identify
# the principal components that explain the most variance. We then
# perform Wilcoxon rank-sum tests to determine whether there are
# significant differences in the top PCs between the two clones,
# and visualise these differences with a scatter plot and a loadings
# heatmap. Plots are saved to project/figures/PCA/.

# ================================================================

# Packages required ----------------------------------------------

# for data manipulation and visualisation
library(tidyverse)

# for interactive plots
library(plotly)

# for reading Excel files
library(readxl)


# Data import ----------------------------------------------------

# filtered livecyte data
metrics <- read_tsv(
  "project/data/movement_morphology/livecyte_collapsed_filtered.tsv",
  show_col_types = FALSE)

# osteogenesis data, pivoted to long format
osteo <- read_xlsx("project/data/differentiation/osteogenesis_processed.xlsx") |>
  pivot_longer(cols = c("0", "8", "8osteo"),
               names_to = "day",
               values_to = "absorption")


# Principal component analysis -----------------------------------

# Select only numeric columns for PCA, removing replicate and
# tracking.id along with n_frames and start_frame which are
# non-biological parameters
numeric_metrics <- metrics |>
  select(-replicate, -tracking.id, -n_frames, -start_frame) |>
  select(where(is.numeric))

# Perform PCA with scaling so all features contribute equally
# regardless of their original units
pca_result <- prcomp(numeric_metrics, scale. = TRUE)

# Create a data frame combining PCA scores with clone labels
pca_df <- as.data.frame(pca_result$x) |>
  bind_cols(metrics |> select(-where(is.numeric)))


# Scree plot -----------------------------------------------------

# Calculate the proportion of variance explained (PVE) by each PC
pve <- pca_result$sdev^2 / sum(pca_result$sdev^2)

data.frame(PC = paste0("PC", 1:length(pve)), PVE = pve) |>
  ggplot(aes(x = PC, y = PVE)) +
  geom_bar(stat = "identity")

# PC1 and PC2 explain ~60% of the variance in the data, so we will
# focus on these two PCs for further analysis. PCs after this become
# less biologically interpretable.


# Wilcoxon rank-sum tests on PCs ---------------------------------

# PC1 is largely loaded with morphological parameters. PC2 also
# considers morphology but additionally loads on total path length
# and final displacement. Neither PC loads heavily on mean speed,
# but mean speed was not significantly different between the clones
# in the LMM analysis, so this is not a concern.

wilcox_PC1 <- wilcox.test(PC1 ~ clone, data = pca_df)
wilcox_PC2 <- wilcox.test(PC2 ~ clone, data = pca_df)

# Both tests are highly significant (p < 0.001), indicating that
# the two clones occupy distinct regions of PC space.

# Extract median values for each clone for annotation
med_PC1_A <- median(pca_df$PC1[pca_df$clone == "A"])
med_PC1_B <- median(pca_df$PC1[pca_df$clone == "B"])
med_PC2_A <- median(pca_df$PC2[pca_df$clone == "A"])
med_PC2_B <- median(pca_df$PC2[pca_df$clone == "B"])


# PC scatter plot ------------------------------------------------

# Build PC2 vs PC1 scatter plot with 95% confidence ellipses and
# median points for each clone annotated with significance brackets
PC2_PC1_plot <- ggplot(pca_df,
                       aes(x = PC1, y = PC2,
                           colour = clone, alpha = clone)) +
  geom_point(size = 3) +
  scale_alpha_manual(values = c("A" = 0.4, "B" = 0.25), guide = "none") +
  labs(x = "PC1", y = "PC2", colour = "Clone") +
  theme_test() +
  stat_ellipse(level = 0.95,
               aes(fill = clone),
               alpha = 0, geom = "polygon", show.legend = FALSE) +
  theme(legend.position = "none") +
  # median point for clone A
  geom_point(aes(x = med_PC1_A, y = med_PC2_A),
             fill = "red", colour = "red", shape = 24, size = 3) +
  # median point for clone B
  geom_point(aes(x = med_PC1_B, y = med_PC2_B),
             colour = "blue", fill = "blue", shape = 24, size = 3) +
  # PC2 significance bracket (vertical, left side)
  annotate("segment",
           x = -4.1, xend = -4.1,
           y = med_PC2_A, yend = med_PC2_B) +
  annotate("segment",
           x = -4.1, xend = -4.1 + 0.1,
           y = med_PC2_A, yend = med_PC2_A) +
  annotate("segment",
           x = -4.1, xend = -4.1 + 0.1,
           y = med_PC2_B, yend = med_PC2_B) +
  annotate("text", x = -4.4, y = -0.2,
           label = "***", size = 4.5) +
  # PC1 significance bracket (horizontal, bottom)
  annotate("segment",
           x = med_PC1_B, xend = med_PC1_A, y = -3.75) +
  annotate("segment",
           x = med_PC1_B, y = -3.75, yend = -3.75 + 0.35) +
  annotate("segment",
           x = med_PC1_A, y = -3.75, yend = -3.75 + 0.35) +
  annotate("text", x = 0, y = -4.2,
           label = "***", size = 4.5)
PC2_PC1_plot


# Loadings heatmap -----------------------------------------------

# Extract loadings for the first 2 PCs to show which features
# contribute most to each component
loadings_mat <- pca_result$rotation[, 1:2]

loadings_long <- as.data.frame(loadings_mat) |>
  rownames_to_column("feature") |>
  pivot_longer(cols = PC1:PC2,
               names_to = "PC",
               values_to = "loading") |>
  mutate(feature = factor(feature, levels = rev(rownames(loadings_mat))))

loadings_heatmap <- ggplot(loadings_long,
                           aes(x = PC, y = feature, fill = loading)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(loading, 2)), size = 3.2) +
  scale_fill_gradient2(low = "#ED7117", mid = "white", high = "#6F2DA8",
                       midpoint = 0, limits = c(-0.6, 0.7)) +
  labs(x = NULL, y = NULL, fill = "Loading") +
  theme_test() +
  theme(legend.position = "right")
loadings_heatmap


# Save plots -----------------------------------------------------

units <- "in"
dpi <- 300
device <- "png"

ggsave("project/figures/PCA/pca_PC2_PC1_plot.png",
       plot = PC2_PC1_plot,
       device = device,
       width = 5, height = 4,
       units = units, dpi = dpi)

ggsave("project/figures/PCA/pca_loadings_heatmap.png",
       plot = loadings_heatmap,
       device = device,
       width = 5, height = 3.75,
       units = units, dpi = dpi)


# R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"
# R Core Team (2024). _R: A Language and Environment for Statistical
# Computing_. R Foundation for Statistical Computing, Vienna, Austria.
# <https://www.R-project.org/>.

# Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R,
# Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL,
# Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP,
# Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019).
# "Welcome to the tidyverse." _Journal of Open Source Software_,
# *4*(43), 1686. doi:10.21105/joss.01686
# <https://doi.org/10.21105/joss.01686>.

# Sievert C (2020). _Interactive Web-Based Data Visualization with R,
# plotly, and shiny_. Chapman and Hall/CRC.
# <https://plotly-r.com>.

# Wickham H, Bryan J (2023). _readxl: Read Excel Files_. R package
# version 1.4.3, <https://CRAN.R-project.org/package=readxl>.