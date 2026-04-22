# ====================================================================
# Livecyte Data Analysis
# - Movement figures:
# - Mean speed violin + box plot
# ====================================================================
#
# 1. Load libraries --------------------------------------------------
#
library(ggplot2) # plotting
library(emmeans) # estimated marginal means (for LMM post-hoc tests)
library(dplyr)  # data manipulation
library(patchwork) # combining plots
library(ggh4x) # for nested facet labels
#
# 2. Load data -------------------------------------------------------
#
df <- read.delim('project/data/movement_morphology/livecyte_collapsed_filtered.tsv')
#
# 3. Define colours and summaries ------------------------------------
#
# We already have p-values from our LMM:
#
# feature               estimate     se ci_lower ci_upper t_value      df   p_value marginal_r2 conditional_r2 replicate_var residual_var     icc p_adjusted sig  
# dry.mass                0.0995 0.0704  -0.0385   0.238     1.41    3.30 0.244         0.00247        0.00489       0.00243        0.997 0.00243  0.279     ""   
# volume                  0.0995 0.0704  -0.0385   0.238     1.41    3.30 0.244         0.00247        0.00489       0.00243        0.997 0.00243  0.279     ""   
# radius                  1.42   0.0547   1.31     1.52     25.9     3.74 0.0000231     0.498          0.500         0.00190        0.503 0.00377  0.0000923 "***"
# sphericity             -1.75   0.0271  -1.80    -1.70    -64.6  1274.   0             0.766          0.766         0              0.234 0        0         "***"
# length.to.width.ratio  -1.31   0.0888  -1.48    -1.13    -14.7     3.98 0.000127      0.421          0.429         0.00866        0.578 0.0148   0.000339  "***"
# mean.speed             -0.437  0.185   -0.800   -0.0730   -2.35    3.85 0.0807        0.0458         0.0904        0.0462         0.944 0.0467   0.129     ""   
# total_path_length      -0.122  0.114   -0.346    0.103    -1.06    4.07 0.347         0.00367        0.0179        0.0142         0.987 0.0142   0.347     ""   
# final_displacement     -0.649  0.108   -0.860   -0.438    -6.02    3.75 0.00469       0.104          0.117         0.0126         0.890 0.0139   0.00937   "**"
# mean.thickness         -1.79   0.0251  -1.84    -1.74    -71.3  1274.   0             0.800          0.800        3.15e-16        0.200 1.57e-15  0        "***"
#
# 3ai. Define clone colour scheme
#
clone_cols   <- c("A" = "#F8766D", "B" = "#00BFC4")
clone_fills  <- c("A" = "#F8766D", "B" = "#00BFC4")
#
# 3aii. Define replicate colour scheme
#
rep_cols <- c(
"A.1" = "#e8564a", "A.2" = "#F8766D", "A.3" = "#f9a090",
"B.1" = "#009ea3", "B.2" = "#00BFC4", "B.3" = "#5dd9dd")
#
# 3b. Combine clone and replicate into a single variable
#
df$clone     <- factor(df$clone, levels = c("A", "B"))
df$replicate <- factor(df$replicate)
df$clone_rep <- interaction(df$clone, df$replicate)
#
# 4. Calculate summaries ---------------------------------------------
#
# 4a. Replicate means
#
rep_means <- df |> 
  group_by(clone, replicate, clone_rep) |> 
  summarise(
    volume = mean(volume,         na.rm = TRUE),
    radius = mean(radius, na.rm = TRUE),
    mean_thickness = mean(mean.thickness,  na.rm = TRUE),
    sphericity = mean(sphericity, na.rm = TRUE),
    length_to_width = mean(length.to.width.ratio, na.rm = TRUE),
    dry_mass = mean(dry.mass, na.rm = TRUE),
    .groups = "drop"
  )
#
# 4b. Clone means + summaries
#
clone_summary <- rep_means |> 
  group_by(clone) |> 
  summarise(
     v_mean = mean(volume), v_se = sd(volume) / sqrt(n()),
     r_mean = mean(radius), r_se = sd(radius) / sqrt(n()),
     t_mean = mean(mean_thickness), t_se = sd(mean_thickness) / sqrt(n()),
     s_mean = mean(sphericity), s_se = sd(sphericity) / sqrt(n()),
     lw_mean = mean(length_to_width), lw_se = sd(length_to_width) / sqrt(n()),
     dm_mean = mean(dry_mass), dm_se = sd(dry_mass) / sqrt(n()),
    .groups = "drop"
  )
#
# 5. Plotting --------------------------------------------------------
#
# 5a. Volume
#
# 5ai. Define label positions
#
v_ymax   <- max(df$volume, na.rm = TRUE)
v_yrange <- v_ymax - min(df$volume, na.rm = TRUE)
v_brack  <- v_ymax + v_yrange * 0.06
v_tick   <- v_yrange * 0.02
v_label  <- v_brack + v_yrange * 0.04
#
# 5aii. Build plot
#
p_volume <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = clone, y = volume, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_crossbar(
    data = clone_summary,
    aes(x = clone, y = v_mean, ymin = v_mean - v_se, ymax = v_mean + v_se),
    width = 0.3, fatten = 1.5, linewidth = 0.4, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = volume, colour = clone_rep),
    size = 3, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2, y = v_brack, yend = v_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = v_brack, yend = v_brack - v_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = v_brack, yend = v_brack - v_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = v_label, label = "p = 0.279", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL, y = expression("Mean Volume (" * mu * "m" ^ 3 * ")")) +
  theme_bw()
p_volume
#
# 5b. Radius
#
# 5bi. Define label positions
#
r_ymax   <- max(df$radius, na.rm = TRUE)
r_yrange <- r_ymax - min(df$radius, na.rm = TRUE)
r_brack  <- r_ymax + r_yrange * 0.06
r_tick   <- r_yrange * 0.02
r_label  <- r_brack + r_yrange * 0.04
#
# 5bii. Build plot
#
p_radius <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = clone, y = radius, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_crossbar(
    data = clone_summary,
    aes(x = clone, y = r_mean, ymin = r_mean - r_se, ymax = r_mean + r_se),
    width = 0.3, fatten = 1.5, linewidth = 0.4, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = radius, colour = clone_rep),
    size = 3, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2, y = r_brack, yend = r_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = r_brack, yend = r_brack - r_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = r_brack, yend = r_brack - r_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = r_label, label = "p = 0.0000923", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL, y = expression("Mean Radius (" * mu * "m)")) +
  theme_bw()
p_radius
#
# 5c. Sphericity
#
# 5ci. Define label positions
#
s_ymax   <- max(df$sphericity, na.rm = TRUE)
s_yrange <- s_ymax - min(df$sphericity, na.rm = TRUE)
s_brack  <- s_ymax + s_yrange * 0.06
s_tick   <- s_yrange * 0.02
s_label  <- s_brack + s_yrange * 0.04
#
# 5cii. Build plot
#
p_sphericity <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = clone, y = sphericity, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_crossbar(
    data = clone_summary,
    aes(x = clone, y = s_mean, ymin = s_mean - s_se, ymax = s_mean + s_se),
    width = 0.3, fatten = 1.5, linewidth = 0.4, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = sphericity, colour = clone_rep),
    size = 3, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2, y = s_brack, yend = s_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = s_brack, yend = s_brack - s_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = s_brack, yend = s_brack - s_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = s_label, label = "p = 0.0000", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL, y = "Sphericity") +
  theme_bw()
p_sphericity
#
# 5d. Length-to-width ratio
#
# 5di. Define label positions
#
lw_ymax   <- max(df$length.to.width.ratio, na.rm = TRUE)
lw_yrange <- lw_ymax - min(df$length.to.width.ratio, na.rm = TRUE)
lw_brack  <- lw_ymax + lw_yrange * 0.06
lw_tick   <- lw_yrange * 0.02
lw_label  <- lw_brack + lw_yrange * 0.04
#
# 5dii. Build plot
#
p_ltwr <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = clone, y = length.to.width.ratio, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_crossbar(
    data = clone_summary,
    aes(x = clone, y = lw_mean, ymin = lw_mean - lw_se, ymax = lw_mean + lw_se),
    width = 0.3, fatten = 1.5, linewidth = 0.4, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = length_to_width, colour = clone_rep),
    size = 3, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2, y = lw_brack, yend = lw_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = lw_brack, yend = lw_brack - lw_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = lw_brack, yend = lw_brack - lw_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = lw_label, label = "p = 0.000339", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL, y = "Length-to-width ratio") +
  theme_bw()
p_ltwr
#
# 5e. Dry mass, mean thickness, and sphericity are all very very highly correlated, as if they're 1:1 transformations, so we will represent this parameter simply as sphericity.
#
# 6. Save all plots --------------------------------------------------
#
ggsave("project/figures/morphology/volume_plot.png", p_volume, width = 4, height = 5, dpi = 300)
ggsave("project/figures/morphology/radius_plot.png", p_radius, width = 4, height = 5, dpi = 300)
ggsave("project/figures/morphology/sphericity_plot.png", p_sphericity, width = 4, height = 5, dpi = 300)
ggsave("project/figures/morphology/length_to_width_plot.png", p_ltwr, width = 4, height = 5, dpi = 300)
