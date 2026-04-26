# Morphology figures for clonal comparison -----------------------

# ================================================================

# Experimental overview ------------------------------------------

# The Livecyte dataset contains measurements of various morphological
# features of cells from two clonal lines (A and B), with three
# replicates each. This script generates the comparison plots for
# volume, radius, sphericity, and length-to-width ratio between the
# two clones.


# Description of data --------------------------------------------

# The filtered livecyte dataset (livecyte_collapsed_filtered.tsv)
# contains one row per tracking.id with columns: clone, replicate,
# tracking.id, and morphological features including volume, radius,
# sphericity, length.to.width.ratio, dry.mass, and mean.thickness.
# The data is stored in project/data/movement_morphology/.


# Analysis overview ----------------------------------------------

# Statistical analyses (LMMs) have already been performed in
# 3.linear_mixed_models.R; the relevant p-values are included here
# for annotation. For each morphological feature, we plot individual
# cell values as jittered points coloured by replicate, replicate
# means as diamonds, and clone-level mean + or - SE as crossbars.
# Significance brackets with BH-adjusted p-values are added.
# Mean thickness, and sphericity are very highly correlated
# (near 1:1 transformations), so we represent this group simply as
# sphericity. Likewise with dry.mass and volume.
# Figures are saved to project/figures/morphology/.

# ================================================================

# Packages required ----------------------------------------------

# for plotting
library(ggplot2)

# for estimated marginal means (LMM post-hoc tests)
library(emmeans)

# for data manipulation
library(dplyr)

# for combining plots
library(patchwork)

# for nested facet labels
library(ggh4x)


# Data import ----------------------------------------------------

df <- read.delim(
  "project/data/movement_morphology/livecyte_collapsed_filtered.tsv")


# LMM p-values for annotation -----------------------------------

# These are BH-adjusted p-values from the LMM analysis in script 3:
#   volume:                p = 0.279
#   radius:                p = 0.0000923  ***
#   sphericity:            p = 0.0000     ***
#   length.to.width.ratio: p = 0.000339   ***
#   dry.mass:              p = 0.279
#   mean.thickness:        p = 0.0000     ***


# Define colours -------------------------------------------------

clone_cols  <- c("A" = "#F8766D", "B" = "#00BFC4")
clone_fills <- c("A" = "#F8766D", "B" = "#00BFC4")

# Replicate-level shades within each clone colour
rep_cols <- c(
  "A.1" = "#e8564a", "A.2" = "#F8766D", "A.3" = "#f9a090",
  "B.1" = "#009ea3", "B.2" = "#00BFC4", "B.3" = "#5dd9dd")


# Prepare data ---------------------------------------------------

df$clone     <- factor(df$clone, levels = c("A", "B"))
df$replicate <- factor(df$replicate)
df$clone_rep <- interaction(df$clone, df$replicate)


# Calculate summaries --------------------------------------------

# Replicate means
rep_means <- df |>
  group_by(clone, replicate, clone_rep) |>
  summarise(
    volume          = mean(volume, na.rm = TRUE),
    radius          = mean(radius, na.rm = TRUE),
    mean_thickness  = mean(mean.thickness, na.rm = TRUE),
    sphericity      = mean(sphericity, na.rm = TRUE),
    length_to_width = mean(length.to.width.ratio, na.rm = TRUE),
    dry_mass        = mean(dry.mass, na.rm = TRUE),
    .groups = "drop"
  )

# Clone means
clone_summary <- df |>
  group_by(clone) |>
  summarise(
    v_mean = mean(volume, na.rm = TRUE),
    r_mean = mean(radius, na.rm = TRUE),
    s_mean = mean(sphericity, na.rm = TRUE),
    lw_mean = mean(length.to.width.ratio, na.rm = TRUE))

# Volume plot ----------------------------------------------------

# Label positions for significance bracket
v_ymax   <- max(df$volume, na.rm = TRUE)
v_yrange <- v_ymax - min(df$volume, na.rm = TRUE)
v_brack  <- v_ymax + v_yrange * 0.06
v_tick   <- v_yrange * 0.02
v_label  <- v_brack + v_yrange * 0.04

p_volume <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = clone, y = volume, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_boxplot(
    data = df,
    aes(x = clone, y = volume),
    width = 0.3, outlier.shape = NA, alpha = 0.5, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = volume, colour = clone_rep),
    size = 3, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2,
           y = v_brack, yend = v_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1,
           y = v_brack, yend = v_brack - v_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2,
           y = v_brack, yend = v_brack - v_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = v_label,
           label = "p = 0.279", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL, y = expression("Volume (" * mu * "m"^3 * ")")) +
  theme_bw()
p_volume


# Radius plot ----------------------------------------------------

r_ymax   <- max(df$radius, na.rm = TRUE)
r_yrange <- r_ymax - min(df$radius, na.rm = TRUE)
r_brack  <- r_ymax + r_yrange * 0.06
r_tick   <- r_yrange * 0.02
r_label  <- r_brack + r_yrange * 0.04

p_radius <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = clone, y = radius, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_boxplot(
    data = df,
    aes(x = clone, y = radius),
    width = 0.3, outlier.shape = NA, alpha = 0.5, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = radius, colour = clone_rep),
    size = 3, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2,
           y = r_brack, yend = r_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1,
           y = r_brack, yend = r_brack - r_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2,
           y = r_brack, yend = r_brack - r_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = r_label,
           label = "p = 0.0000923", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL, y = expression("Radius (" * mu * "m)")) +
  theme_bw()
p_radius


# Sphericity plot ------------------------------------------------

s_ymax   <- max(df$sphericity, na.rm = TRUE)
s_yrange <- s_ymax - min(df$sphericity, na.rm = TRUE)
s_brack  <- s_ymax + s_yrange * 0.06
s_tick   <- s_yrange * 0.02
s_label  <- s_brack + s_yrange * 0.04

p_sphericity <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = clone, y = sphericity, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_boxplot(
    data = df,
    aes(x = clone, y = sphericity),
    width = 0.3, outlier.shape = NA, alpha = 0.5, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = sphericity, colour = clone_rep),
    size = 3, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2,
           y = s_brack, yend = s_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1,
           y = s_brack, yend = s_brack - s_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2,
           y = s_brack, yend = s_brack - s_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = s_label,
           label = "p < 0.00001", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL, y = "Sphericity") +
  theme_bw()
p_sphericity


# Length-to-width ratio plot -------------------------------------

lw_ymax   <- max(df$length.to.width.ratio, na.rm = TRUE)
lw_yrange <- lw_ymax - min(df$length.to.width.ratio, na.rm = TRUE)
lw_brack  <- lw_ymax + lw_yrange * 0.06
lw_tick   <- lw_yrange * 0.02
lw_label  <- lw_brack + lw_yrange * 0.04

p_ltwr <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = clone, y = length.to.width.ratio, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_boxplot(
    data = df,
    aes(x = clone, y = length.to.width.ratio),
    width = 0.3, outlier.shape = NA, alpha = 0.5, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = length_to_width, colour = clone_rep),
    size = 3, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2,
           y = lw_brack, yend = lw_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1,
           y = lw_brack, yend = lw_brack - lw_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2,
           y = lw_brack, yend = lw_brack - lw_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = lw_label,
           label = "p = 0.000339", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL, y = "Aspect ratio") +
  theme_bw()
p_ltwr

# Dry mass, mean thickness, and sphericity are all very highly
# correlated, as if they're 1:1 transformations, so we'll represent
# this parameter simply as sphericity.


# Save all plots -------------------------------------------------

ggsave("project/figures/morphology/volume_plot.png",
       plot = p_volume, width = 4, height = 5, dpi = 300)

ggsave("project/figures/morphology/radius_plot.png",
       plot = p_radius, width = 4, height = 5, dpi = 300)

ggsave("project/figures/morphology/sphericity_plot.png",
       plot = p_sphericity, width = 4, height = 5, dpi = 300)

ggsave("project/figures/morphology/length_to_width_plot.png",
       plot = p_ltwr, width = 4, height = 5, dpi = 300)


# R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"
# R Core Team (2024). _R: A Language and Environment for Statistical
# Computing_. R Foundation for Statistical Computing, Vienna, Austria.
# <https://www.R-project.org/>.

# Wickham H (2016). _ggplot2: Elegant Graphics for Data Analysis_.
# Springer-Verlag New York. <https://ggplot2.tidyverse.org>.

# Lenth R (2024). _emmeans: Estimated Marginal Means, aka
# Least-Squares Means_. R package version 1.10.5,
# <https://CRAN.R-project.org/package=emmeans>.

# Wickham H, François R, Henry L, Müller K, Vaughan D (2023).
# _dplyr: A Grammar of Data Manipulation_. R package version 1.1.4,
# <https://CRAN.R-project.org/package=dplyr>.

# Pedersen TL (2024). _patchwork: The Composer of Plots_. R package
# version 1.3.0, <https://CRAN.R-project.org/package=patchwork>.

# van den Brand T (2024). _ggh4x: Hacks for 'ggplot2'_. R package
# version 0.2.8, <https://CRAN.R-project.org/package=ggh4x>.