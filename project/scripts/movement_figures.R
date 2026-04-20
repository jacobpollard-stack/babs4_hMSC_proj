# ==========================================================================
# Livecyte Data Analysis
# - Movement figures:
# - Mean speed violin + box plot
# ==========================================================================
#
# 1. Load libraries
#
library(lme4)
library(lmerTest)
library(ggplot2)
library(emmeans)
library(dplyr)
library(patchwork)
library(ggh4x)
#
# 2. Load data
#
df <- read.delim("project/data/movement_morphology/livecyte_collapsed_filtered.tsv")
#
collapsed <- read.delim('project/data/movement_morphology/livecyte_filtered.tsv')
#
# 3. Define colours and summaries
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
#
# 3ai. Define colour scheme
#
clone_cols   <- c("A" = "#F8766D", "B" = "#00BFC4")
clone_fills  <- c("A" = "#F8766D", "B" = "#00BFC4")
#
# 3aii. Define replicate shades
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
# 3c. Calculate replicate means
#
rep_means <- df |> 
  group_by(clone, replicate, clone_rep) |> 
  summarise(
    mean_speed = mean(mean.speed,         na.rm = TRUE),
    final_disp = mean(final_displacement, na.rm = TRUE),
    total_path = mean(total_path_length,  na.rm = TRUE),
    .groups = "drop"
  )
#
# 3d. Calculate clone means + summaries
#
clone_summary <- rep_means |> 
  group_by(clone) |> 
  summarise(
    ms_mean = mean(mean_speed), ms_se = sd(mean_speed) / sqrt(n()),
    fd_mean = mean(final_disp), fd_se = sd(final_disp) / sqrt(n()),
    tp_mean = mean(total_path), tp_se = sd(total_path) / sqrt(n()),
    .groups = "drop"
  )
#
# 4. Form mean speed plot
#
# 4a. Define label positions
#
ms_ymax   <- max(df$mean.speed, na.rm = TRUE)
ms_yrange <- ms_ymax - min(df$mean.speed, na.rm = TRUE)
ms_brack  <- ms_ymax + ms_yrange * 0.06
ms_tick   <- ms_yrange * 0.02
ms_label  <- ms_brack + ms_yrange * 0.04
#
# 4b. Mean speed plot
#
p_mean.speed <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = clone, y = mean.speed, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_crossbar(
    data = clone_summary,
    aes(x = clone, y = ms_mean, ymin = ms_mean - ms_se, ymax = ms_mean + ms_se),
    width = 0.3, fatten = 1.5, linewidth = 0.4, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = mean_speed, colour = clone_rep),
    size = 3, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2, y = ms_brack, yend = ms_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = ms_brack, yend = ms_brack - ms_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = ms_brack, yend = ms_brack - ms_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = ms_label, label = "p = 0.081", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL, y = expression("Mean speed (" * mu * "m s"^{-1} * ")")) +
  theme_bw()
p_mean.speed
#
# 5. Final displacement plot
#
# 5a. Define label positions
fd_ymax   <- max(df$final_displacement, na.rm = TRUE)
fd_yrange <- fd_ymax - min(df$final_displacement, na.rm = TRUE)
fd_brack  <- fd_ymax + fd_yrange * 0.06
fd_tick   <- fd_yrange * 0.02
fd_label  <- fd_brack + fd_yrange * 0.04
#
# 5b. Final displacement plot
p_final_displacement <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = clone, y = final_displacement, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_crossbar(
    data = clone_summary,
    aes(x = clone, y = fd_mean, ymin = fd_mean - fd_se, ymax = fd_mean + fd_se),
    width = 0.3, fatten = 1.5, linewidth = 0.4, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = final_disp, colour = clone_rep),
    size = 3.5, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2, y = fd_brack, yend = fd_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = fd_brack, yend = fd_brack - fd_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = fd_brack, yend = fd_brack - fd_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = fd_label, label = "p = 0.00469", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 15) +
  labs(x = NULL, y = expression("Final displacement (" * mu * "m)")) +
  theme_bw()
p_final_displacement
#
# 6. Total path length plot
#
# 6a. Define label positions
#
tp_ymax   <- max(df$total_path_length, na.rm = TRUE)
tp_yrange <- tp_ymax - min(df$total_path_length, na.rm = TRUE)
tp_brack  <- tp_ymax + tp_yrange * 0.06
tp_tick   <- tp_yrange * 0.02
tp_label  <- tp_brack + tp_yrange * 0.04
#
# 6b. Total path length plot
#
p_total_path_length <- ggplot() +
  geom_jitter(
    data = df,
    aes(x = clone, y = total_path_length, colour = clone_rep),
    width = 0.25, size = 0.6, alpha = 0.35, shape = 16
  ) +
  geom_crossbar(
    data = clone_summary,
    aes(x = clone, y = tp_mean, ymin = tp_mean - tp_se, ymax = tp_mean + tp_se),
    width = 0.3, fatten = 1.5, linewidth = 0.4, colour = "black", fill = NA
  ) +
  geom_point(
    data = rep_means,
    aes(x = clone, y = total_path, colour = clone_rep),
    size = 3.5, shape = 18,
    position = position_dodge(width = 0.2)
  ) +
  annotate("segment", x = 1, xend = 2, y = tp_brack, yend = tp_brack, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = tp_brack, yend = tp_brack - tp_tick, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = tp_brack, yend = tp_brack - tp_tick, linewidth = 0.4) +
  annotate("text", x = 1.5, y = tp_label, label = "p = 0.347", size = 3.5) +
  scale_colour_manual(values = rep_cols, guide = "none") +
  scale_x_discrete(labels = c("Clone A", "Clone B")) +
  scale_y_continuous(n.breaks = 10) +
  labs(x = NULL, y = expression("Total path length (" * mu * "m)")) +
  theme_bw()
p_total_path_length
#
# 7. Cell movement plot
#
# 7a. Normalise position.x and position.y values
#
collapsed <- collapsed |> 
  group_by(clone, replicate, tracking.id) |> 
  arrange(frame, .by_group = TRUE) |> 
  mutate(
    dx = position.x - first(position.x),
    dy = position.y - first(position.y)
  ) |> 
  ungroup()
#
# 7b. Form track plot
#
# 7bi. Sample size per clone = 50
#
set.seed(42)
#
sampled_ids <- collapsed |> 
  distinct(clone, replicate, tracking.id) |> 
  group_by(clone) |> 
  slice_sample(n = 50) |> 
  ungroup()
#
df_sample <- collapsed |> 
  semi_join(sampled_ids, by = c("clone", "replicate", "tracking.id"))
#
# 7bii. Rename A and B to clone A and clone B for facet header
#
df_sample$clone <- factor(df_sample$clone,
                    levels = c("A", "B"),
                    labels = c("Clone A", "Clone B"))
#
# 7biii. Plot
#
endpoints <- df_sample |> 
  group_by(clone, replicate, tracking.id) |> 
  slice_tail(n = 1) |> 
  ungroup()
#
p_spaghetti <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey", linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey", linewidth = 0.4, linetype = "dashed") +
  geom_path(
    data = df_sample,
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
#
# 8. Save all plots
#
ggsave("project/figures/movement/mean.speed_plot.png", p_mean.speed, width = 4, height = 5, dpi = 300)
ggsave("project/figures/movement/final_displacement_plot.png", p_final_displacement, width = 4, height = 5, dpi = 300)
ggsave("project/figures/movement/total_path_length_plot.png", p_total_path_length, width = 4, height = 5, dpi = 300)
ggsave("project/figures/movement/tracking_figure.png", p_spaghetti, width = 4, height = 5, dpi = 300)

