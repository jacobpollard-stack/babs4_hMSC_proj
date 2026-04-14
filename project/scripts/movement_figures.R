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
# 3. Statistical analysis of movement data
#
# 3a. We have hierarchal data, so we will use a Linear Mixed Model
# We have 2 clones, and 3 replicates per clone, so we will use clone as a fixed effect and replicate as a random effect.
# We will follow the LMM with ANOVA to test for significance of the fixed effect.
#
# 3ai. mean.speed
#
lmm_mean.speed <- lmer(mean.speed ~ clone + (1 | clone:replicate), data = df)
mean.speed_anova <- anova(lmm_mean.speed)
print(mean.speed_anova)
#
# Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq  Mean Sq NumDF  DenDF F value  Pr(>F)  
# clone 0.078854 0.078854     1 3.8531  5.5403 0.08068 .
#
# There is a trend towards significance for the effect of clone on mean speed, but it does not reach significance at the 0.05 level (F(1, 3.85) = 5.54, p = 0.081).
#
# 3aii. final_displacement
#
lmm_final_displacement <- lmer(final_displacement ~ clone + (1 | clone:replicate), data = df)
final_displacement_anova <- anova(lmm_final_displacement)
#
# Type III Analysis of Variance Table with Satterthwaite's method
#        Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# clone 163150  163150     1 3.7544  36.252 0.004687 **
#
# There is a significant effect of clone on final displacement (F(1, 3.75) = 36.25, p = 0.004687), indicating that the two clones may differ in their final displacement: clone A had a higher median final displacement (129.6μm) compared to clone B (91.0μm).
#
# 3aiii. total_path_length
#
lmm_total_path_length <- lmer(total_path_length ~ clone + (1 | clone:replicate), data = df)
total_path_length_anova <- anova(lmm_total_path_length)
#
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
# clone  20874   20874     1 4.0714  1.1277 0.3471
#
# There is no significant effect of clone on total path length (F(1, 4.07) = 1.13, p = 0.3471), indicating that the two clones do not differ in their total path length.
#
# 4. Plotting movement data
#
# 4ai. Define colour scheme
#
clone_cols   <- c("A" = "#F8766D", "B" = "#00BFC4")
clone_fills  <- c("A" = "#F8766D", "B" = "#00BFC4")
#
# 4aii. Define replicate shades
#
rep_cols <- c(
  "A.1" = "#e8564a", "A.2" = "#F8766D", "A.3" = "#f9a090",
  "B.1" = "#009ea3", "B.2" = "#00BFC4", "B.3" = "#5dd9dd")
#
# 4b. Combine clone and replicate into a single variable
#
df$clone     <- factor(df$clone, levels = c("A", "B"))
df$replicate <- factor(df$replicate)
df$clone_rep <- interaction(df$clone, df$replicate)
#
# 4c. Calculate replicate means
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
# 4d. Calculate clone means + summaries
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
# 5. Form mean speed plot
#
# 5a. Define label positions
#
ms_ymax   <- max(df$mean.speed, na.rm = TRUE)
ms_yrange <- ms_ymax - min(df$mean.speed, na.rm = TRUE)
ms_brack  <- ms_ymax + ms_yrange * 0.06
ms_tick   <- ms_yrange * 0.02
ms_label  <- ms_brack + ms_yrange * 0.04
#
# 5b. Mean speed plot
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
# 6. Final displacement plot
#
# 6a. Define label positions
fd_ymax   <- max(df$final_displacement, na.rm = TRUE)
fd_yrange <- fd_ymax - min(df$final_displacement, na.rm = TRUE)
fd_brack  <- fd_ymax + fd_yrange * 0.06
fd_tick   <- fd_yrange * 0.02
fd_label  <- fd_brack + fd_yrange * 0.04
#
# 6b. Final displacement plot
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
# 7. Total path length plot
#
# 7a. Define label positions
#
tp_ymax   <- max(df$total_path_length, na.rm = TRUE)
tp_yrange <- tp_ymax - min(df$total_path_length, na.rm = TRUE)
tp_brack  <- tp_ymax + tp_yrange * 0.06
tp_tick   <- tp_yrange * 0.02
tp_label  <- tp_brack + tp_yrange * 0.04
#
# 7b. Total path length plot
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
# 8. Cell movement plot
#
# 8a. Normalise position.x and position.y values
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
# 8b. Form track plot
#
# 8bi. Sample size per clone = 50
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
# 8bii. Rename A and B to clone A and clone B for facet header
#
df_sample$clone <- factor(df_sample$clone,
                    levels = c("A", "B"),
                    labels = c("Clone A", "Clone B"))
#
# 8biii. Plot
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
# 9. Save all plots
#
ggsave("project/figures/movement_morphology/mean.speed_plot.png", p_mean.speed, width = 4, height = 5, dpi = 300)
ggsave("project/figures/movement_morphology/final_displacement_plot.png", p_final_displacement, width = 4, height = 5, dpi = 300)
ggsave("project/figures/movement_morphology/total_path_length_plot.png", p_total_path_length, width = 4, height = 5, dpi = 300)
ggsave("project/figures/movement_morphology/tracking_figure.png", p_spaghetti, width = 4, height = 5, dpi = 300)

