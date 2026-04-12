# ==========================================================================
# Osteogenesis Differentiation Data Analysis
# ==========================================================================
#
# 1. Load libraries --------------------------------------------------------
#
library(tidyverse)
library(rcompanion)
library(dunn.test)
library(ggh4x)
#
# 2. Load data -------------------------------------------------------------
#
osteo  <- read_xlsx("project/data/differentiation/osteogenesis_processed.xlsx")
#
# 3. Pivot longer ----------------------------------------------------------
#
osteo_long <- osteo |> 
  pivot_longer(cols = c("0", "8", "8osteo"), names_to = "day", values_to = "absorption") |> 
  mutate(day = as.factor(day))
#
# 4. Statistical tests -----------------------------------------------------
#
# 4a. Clone A vs clone B at day 0
#
# 4ai. Subset data for day 0 to quantify the effect of clone at baseline
#
osteo_long_day0 <- osteo_long |>
  filter(day == "0")
#
# 4aii. Perform wilcoxon rank-sum test
#
wilcox_AvB_day0 <- wilcox.test(absorption ~ clone, data = osteo_long_day0, exact = FALSE) # as some data points are identical
#
# The Wilcoxon rank-sum test results show no significant difference in absorption between clone A and clone B at day 0 (p = 0.0765225), suggesting that the two clones have similar baseline absorption levels before the introduction of osteogenic medium. However it seems as this is due to low sample size (n = 3). |Δ| ≈ 0.05, which is quite large for this data. Thus it seems that the baseline absorption is different between the clones, but the statistical power to detect this difference is very low.
#
# 4b. Day 8 vs day 0 within each clone
#
# 4bi. Clone A
#
wilcox_A_day8v0 <- wilcox.test(
  osteo_long$absorption[osteo_long$clone == "A" & osteo_long$day == "0"],
  osteo_long$absorption[osteo_long$clone == "A" & osteo_long$day == "8"], exact = FALSE)
#
# 4bii. Clone B
#
wilcox_B_day8v0 <- wilcox.test(
  osteo_long$absorption[osteo_long$clone == "B" & osteo_long$day == "0"],
  osteo_long$absorption[osteo_long$clone == "B" & osteo_long$day == "8"], exact = FALSE)
#
# The Wilcoxon rank-sum test results show that the differences in absorption between day 0 and day 8 for both clones are not statistically significant (clone A: p = 0.0765225; clone B: p = 0.1). However this is also likely due to small sample size, as the data does show an increase in absorption (Δ(A) ≈ 0.17, Δ(B) ≈ 0.2).
#
# 4c. Wilcoxon test on the change in absorption between day 8 and day 0 between clones
#
# 4ci. Calculate the change in absorption between day 8 and day 0 for each clone
#
day_change <- osteo_long |> 
  group_by(clone, replicate) |> 
  summarise(change = absorption[day == "8"] - absorption[day == "0"])
#
# 4cii. Perform Wilcoxon test on the change in absorption
#
wilcox_AvB_0to8 <- wilcox.test(change ~ clone, data = day_change)
#
# The Wilcoxon rank-sum test results show no significant difference in the change in absorption between clone A and clone B (p = 0.1), suggesting that the overall change in ALP activity from day 0 to day 8 is not significantly different between the two clones. Δ for B looks slightly higher than Δ for A, but we will assume statistical non-significance for ease of the next step.
#
# 4d. Test change in absorption between day 8 with osteogenic medium and day 8 control between clones
#
# 4di. Calculate the change in absorption between day 8 with osteogenic medium and day 8 control for each clone
#
osteo_change <- osteo_long |>
  group_by(clone, replicate) |> 
  summarise(change = absorption[day == "8osteo"] - absorption[day == "8"])
#
# 4dii. Perform Wilcoxon test on the change in absorption between clone A and clone B
#
wilcox_osteo_change <- wilcox.test(change ~ clone, data = osteo_change)
#
# The Wilcoxon rank-sun test shows no significant difference between the change in the two clone' absorption values between day 8 with osteogenic medium and day 8 control (p = 0.1). However, looking at the data, the difference between the clones is about 0.196, which is quite large.
#
# 5. Plot the osteogenesis data
#
osteo_plot <- ggplot(osteo_long, aes(x = day, y = absorption, group = replicate)) +
  geom_line(aes(colour = clone), alpha = 0.4, linewidth = 0.7) +
  geom_point(aes(colour = clone), size = 2.5) +
  facet_wrap(~ clone) +
  scale_colour_manual(values = c(
    "A" = "#f7766f",
    "B" = "#0dc1c5"),
    guide = "none") +
  labs(
    x = "Day",
    y = "Absorption") +
  theme(legend.position = "none",
    strip.text = element_text(face = "bold")) +
  scale_x_discrete(labels = c(
    "0" = "0",
    "8" = "8",
    "8osteo" = "8 in osteogenic\nmedium")) +
  scale_y_continuous(n.breaks = 10) +
  facet_wrap2(~ clone, strip = strip_themed(
    background_x = list(
      element_rect(fill = "#f7766f"),
      element_rect(fill = "#0dc1c5")),
    text_x = list(
      element_text(colour = "black", face = "bold"),
      element_text(colour = "black", face = "bold")))) +
  theme_bw()
#
# 6. Save plot
#
ggsave("project/figures/osteogenesis/osteo_plot.png", osteo_plot)

