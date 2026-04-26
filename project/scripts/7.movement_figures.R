# Movement figures for clonal comparison -------------------------

# ================================================================

# Experimental overview ------------------------------------------

# The Livecyte dataset contains movement-related measurements for
# cells from two clonal lines (A and B), with three replicates each.
# This script generates figures for mean speed, final displacement,
# total path length, cell tracking trajectories, and mean squared
# displacement (MSD) analysis.


# Description of data --------------------------------------------

# Two filtered datasets are used:
#    livecyte_collapsed_filtered.tsv — one row per tracking.id with
#       summary metrics including mean.speed, final_displacement,
#       and total_path_length
#    livecyte_data_filtered.tsv — per-frame data with position.x,
#       position.y, and frame for trajectory plotting
# Both are stored in project/data/movement_morphology/.


# Analysis overview ----------------------------------------------

# Statistical analyses (LMMs) have already been performed in
# 3.linear_mixed_models.R; the relevant p-values are included here
# for annotation. For each movement feature, we plot individual cell
# values as jittered points coloured by replicate, replicate means as
# diamonds, and clone-level mean + or - SE as crossbars. We also create
# a spaghetti plot of 50 randomly sampled cell trajectories per clone,
# and a log-log MSD plot with alpha exponent fits to characterise the
# nature of cell movement. Figures are saved to
# project/figures/movement/.

# ================================================================

# Packages required ----------------------------------------------

# for plotting
library(ggplot2)

# for estimated marginal means
library(emmeans)

# for data manipulation
library(dplyr)

# for combining plots
library(patchwork)

# for strip_themed in spaghetti plot
library(ggh4x)


# Data import ----------------------------------------------------

livecyte_collapsed_filtered <- read.delim(
  "project/data/movement_morphology/livecyte_collapsed_filtered.tsv")

livecyte_data_filtered <- read.delim(
  "project/data/movement_morphology/livecyte_data_filtered.tsv")


# LMM p-values for annotation -----------------------------------

# These are BH-adjusted p-values from the LMM analysis in script 3:
#   mean.speed:          p = 0.129
#   total_path_length:   p = 0.347
#   final_displacement:  p = 0.00937  **


# Define colours -------------------------------------------------

clone_cols  <- c("A" = "#F8766D", "B" = "#00BFC4")
clone_fills <- c("A" = "#F8766D", "B" = "#00BFC4")

rep_cols <- c(
  "A.1" = "#e8564a", "A.2" = "#F8766D", "A.3" = "#f9a090",
  "B.1" = "#009ea3", "B.2" = "#00BFC4", "B.3" = "#5dd9dd")


# Prepare data ---------------------------------------------------

livecyte_collapsed_filtered$clone     <- factor(
  livecyte_collapsed_filtered$clone, levels = c("A", "B"))
livecyte_collapsed_filtered$replicate <- factor(
  livecyte_collapsed_filtered$replicate)
livecyte_collapsed_filtered$clone_rep <- interaction(
  livecyte_collapsed_filtered$clone,
  livecyte_collapsed_filtered$replicate)


# Calculate summaries --------------------------------------------

# Replicate means
rep_means <- livecyte_collapsed_filtered |>
  group_by(clone, replicate, clone_rep) |>
  summarise(
    mean_speed = mean(mean.speed, na.rm = TRUE),
    final_disp = mean(final_displacement, na.rm = TRUE),
    total_path = mean(total_path_length, na.rm = TRUE),
    .groups = "drop"
  )

# Clone means
clone_summary <- rep_means |>
  group_by(clone) |>
  summarise(
    ms_mean = mean(mean_speed),
    fd_mean = mean(final_disp),
    tp_mean = mean(total_path),
    .groups = "drop"
  )


# Mean speed plot ------------------------------------------------

ms_ymax   <- max(livecyte_collapsed_filtered$mean.speed, na.rm = TRUE)
ms_yrange <- ms_ymax - min(livecyte_collapsed_filtered$mean.speed, na.rm = TRUE)
ms_brack  <- ms_ymax + ms_yrange * 0.06
ms_tick   <- ms_yrange * 0.02
ms_label  <- ms_brack + ms_yrange * 0.04

p_mean.speed <- ggplot() +
  geom_jitter(
    data = livecyte_collapsed_filtered,
    aes(x = clone, y = mean.speed, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_boxplot(
    data = livecyte_collapsed_filtered,
    aes(x = clone, y = mean.speed),
    width = 0.3, outlier.shape = NA, alpha = 0.5, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = mean_speed, colour = clone_rep),
    size = 3, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2,
           y = ms_brack, yend = ms_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1,
           y = ms_brack, yend = ms_brack - ms_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2,
           y = ms_brack, yend = ms_brack - ms_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = ms_label,
           label = "p = 0.081", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL,
       y = expression("Mean speed (" * mu * "m s"^{-1} * ")")) +
  theme_bw()
p_mean.speed


# Final displacement plot ----------------------------------------

fd_ymax   <- max(livecyte_collapsed_filtered$final_displacement, na.rm = TRUE)
fd_yrange <- fd_ymax - min(livecyte_collapsed_filtered$final_displacement, na.rm = TRUE)
fd_brack  <- fd_ymax + fd_yrange * 0.06
fd_tick   <- fd_yrange * 0.02
fd_label  <- fd_brack + fd_yrange * 0.04

p_final_displacement <- ggplot() +
  geom_jitter(
    data = livecyte_collapsed_filtered,
    aes(x = clone, y = final_displacement, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_boxplot(
    data = livecyte_collapsed_filtered,
    aes(x = clone, y = final_displacement),
    width = 0.3, outlier.shape = NA, alpha = 0.5, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = final_disp, colour = clone_rep),
    size = 3.5, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2,
           y = fd_brack, yend = fd_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1,
           y = fd_brack, yend = fd_brack - fd_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2,
           y = fd_brack, yend = fd_brack - fd_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = fd_label,
           label = "p = 0.00469", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL,
       y = expression("Final displacement (" * mu * "m)")) +
  theme_bw()
p_final_displacement


# Total path length plot -----------------------------------------

tp_ymax   <- max(livecyte_collapsed_filtered$total_path_length, na.rm = TRUE)
tp_yrange <- tp_ymax - min(livecyte_collapsed_filtered$total_path_length, na.rm = TRUE)
tp_brack  <- tp_ymax + tp_yrange * 0.06
tp_tick   <- tp_yrange * 0.02
tp_label  <- tp_brack + tp_yrange * 0.04

p_total_path_length <- ggplot() +
  geom_jitter(
    data = livecyte_collapsed_filtered,
    aes(x = clone, y = total_path_length, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_boxplot(
    data = livecyte_collapsed_filtered,
    aes(x = clone, y = total_path_length),
    width = 0.3, outlier.shape = NA, alpha = 0.5, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = total_path, colour = clone_rep),
    size = 3.5, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2,
           y = tp_brack, yend = tp_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1,
           y = tp_brack, yend = tp_brack - tp_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2,
           y = tp_brack, yend = tp_brack - tp_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = tp_label,
           label = "p = 0.347", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 10) +
  labs(x = NULL,
       y = expression("Total path length (" * mu * "m)")) +
  theme_bw()
p_total_path_length


# Cell tracking spaghetti plot -----------------------------------

# Normalise position.x and position.y so all tracks start at (0, 0)
livecyte_data_filtered <- livecyte_data_filtered |>
  group_by(clone, replicate, tracking.id) |>
  arrange(.data$frame, .by_group = TRUE) |>
  mutate(
    dx = position.x - first(position.x),
    dy = position.y - first(position.y)
  ) |>
  ungroup()

# Randomly sample 50 cells per clone for readability
set.seed(42)

sampled_ids <- livecyte_data_filtered |>
  distinct(clone, replicate, tracking.id) |>
  group_by(clone) |>
  slice_sample(n = 50) |>
  ungroup()

livecyte_data_filtered_sample <- livecyte_data_filtered |>
  semi_join(sampled_ids, by = c("clone", "replicate", "tracking.id"))

# Relabel clones for facet header
livecyte_data_filtered_sample$clone <- factor(
  livecyte_data_filtered_sample$clone,
  levels = c("A", "B"),
  labels = c("Clone A", "Clone B"))

# Extract track endpoints for plotting as dots
endpoints <- livecyte_data_filtered_sample |>
  group_by(clone, replicate, tracking.id) |>
  slice_tail(n = 1) |>
  ungroup()

p_spaghetti <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey",
             linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey",
             linewidth = 0.4, linetype = "dashed") +
  geom_path(
    data = livecyte_data_filtered_sample,
    aes(x = dx, y = dy,
        group = interaction(replicate, tracking.id),
        colour = clone),
    alpha = 0.35, linewidth = 0.4
  ) +
  geom_point(
    data = endpoints,
    aes(x = dx, y = dy, colour = clone),
    size = 0.8, alpha = 0.5
  ) +
  annotate("point", x = 0, y = 0, shape = 3, size = 2, colour = "black") +
  facet_wrap2(~ clone, strip = strip_themed(
    background_x = list(
      element_rect(fill = "#f7766f"),
      element_rect(fill = "#0dc1c5")),
    text_x = list(
      element_text(colour = "black", face = "bold"),
      element_text(colour = "black", face = "bold")))) +
  scale_colour_manual(values = clone_cols, guide = "none") +
  coord_equal() +
  labs(
    x = expression("x displacement (" * mu * "m)"),
    y = expression("y displacement (" * mu * "m)")
  ) +
  theme_classic()
p_spaghetti


# Mean squared displacement (MSD) analysis ----------------------

# Compute MSD for each cell at each time lag. The frame interval is
# 23 minutes (23/60 hours). We cap the maximum lag at 96 hours
frame_interval <- 23 / 60
max_lag_hours <- 96
max_lag <- floor(max_lag_hours / (4 * frame_interval))

msd_cell <- livecyte_data_filtered |>
  group_by(clone, replicate, tracking.id) |>
  arrange(frame, .by_group = TRUE) |>
  do({
    pos <- .
    n   <- nrow(pos)
    lags <- seq_len(min(max_lag, n - 1))
    tibble(
      lag = lags,
      msd = sapply(lags, function(tau) {
        dx <- pos$position.x[(1 + tau):n] - pos$position.x[1:(n - tau)]
        dy <- pos$position.y[(1 + tau):n] - pos$position.y[1:(n - tau)]
        mean(dx^2 + dy^2)
      })
    )
  }) |>
  ungroup()

# Convert lag from frames to hours
msd_cell$time <- msd_cell$lag * frame_interval

# Summarise: replicate-level means, then clone-level mean +/- SE
msd_rep <- msd_cell |>
  group_by(clone, replicate, time) |>
  summarise(msd_mean = mean(msd, na.rm = TRUE), .groups = "drop")

msd_clone <- msd_rep |>
  group_by(clone, time) |>
  summarise(
    msd_grand = mean(msd_mean),
    msd_se    = sd(msd_mean) / sqrt(n()),
    .groups   = "drop"
  )

# MSD plot on log-log scale. Log scale is chosen so the low-lag data
# is more visible, which is important for estimating alpha as it is
# more reliable due to having more data points
p_msd <- ggplot(msd_clone,
                aes(x = time, y = msd_grand,
                    colour = clone, fill = clone)) +
  geom_ribbon(
    aes(ymin = msd_grand - msd_se, ymax = msd_grand + msd_se),
    alpha = 0.2, colour = NA
  ) +
  geom_line(linewidth = 0.8) +
  labs(
    x      = expression(tau * " (h)"),
    y      = expression("MSD (" * mu * "m"^2 * ")"),
    colour = NULL
  ) +
  scale_x_log10(breaks = c(1, 10)) +
  scale_y_log10() +
  scale_colour_manual(values = clone_cols,
                      labels = c("Clone A", "Clone B")) +
  scale_fill_manual(values = clone_fills, guide = "none") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.85)) +
  annotation_logticks(sides = "bl")
p_msd


# Estimate alpha exponent ----------------------------------------

# Fit linear model to log-log MSD data to estimate the anomalous
# diffusion exponent alpha. alpha > 1 indicates superdiffusion,
# alpha = 1 is normal diffusion, alpha < 1 is subdiffusion
alpha_fits <- msd_clone |>
  group_by(clone) |>
  summarise(
    fit       = list(lm(log10(msd_grand) ~ log10(time))),
    alpha     = coef(fit[[1]])[2],
    intercept = coef(fit[[1]])[1],
    .groups = "drop"
  )
alpha_fits


# clone fit    alpha intercept
#  1 A  <lm>    1.28      2.92
#  2 B  <lm>    1.17      2.73
#
# The alpha values suggest that both clones exhibit superdiffusive
# behaviour (alpha > 1), with Clone A being slightly more
# superdiffusive than Clone B.


# Save all plots -------------------------------------------------

ggsave("project/figures/movement/mean.speed_plot.png",
       plot = p_mean.speed,
       width = 4, height = 5, dpi = 300)

ggsave("project/figures/movement/final_displacement_plot.png",
       plot = p_final_displacement,
       width = 4, height = 5, dpi = 300)

ggsave("project/figures/movement/total_path_length_plot.png",
       plot = p_total_path_length,
       width = 4, height = 5, dpi = 300)

ggsave("project/figures/movement/tracking_plot.png",
       plot = p_spaghetti,
       width = 4, height = 5, dpi = 300)

ggsave("project/figures/movement/msd_plot.png",
       plot = p_msd,
       width = 6, height = 4, dpi = 300)


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