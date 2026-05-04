# Osteogenesis differentiation statistical analysis --------------

# ================================================================

# Experimental overview ------------------------------------------

# Three replicates per clonal line were incubated with osteogenic
# differentiation medium for 8 days, and ALP activity was measured at
# day 0 and day 8 with osteogenic medium. DNA concentration was also
# measured in cell lysates at the same time points to control for cell
# number. ALP activity was measured by proxy as absorbance at 405nm.


# Description of data --------------------------------------------

# The osteogenesis data (osteogenesis_processed.xlsx) contains
# absorbance values for each clone and replicate at three time points:
# day 0, day 8, and day 8 with osteogenic medium. The DNA
# concentration data (DNA_cell_lysates.xlsx) contains DNA
# concentration values for each clone and replicate. Both files are
# stored in project/data/differentiation/.


# Analysis overview ----------------------------------------------

# This script performs Wilcoxon rank-sum tests to compare ALP activity
# between the two clones at different time points, and between time
# points within each clone. We also test whether the change in ALP
# activity from day 0 to day 8, and from day 8 control to day 8 with
# osteogenic medium, differs between clones. The results are
# visualised in a plot with absorbance and DNA concentration on dual
# axes. The figure is saved to project/figures/osteogenesis/.

# ================================================================

# Packages required ----------------------------------------------

# for data manipulation and visualiaation
library(tidyverse)

# for the wilcoxon test
library(rcompanion)

# for post-hoc Dunn's test
library(dunn.test)

# for facet_wrap2 with themed strips
library(ggh4x)


# Data import ----------------------------------------------------

# osteogenesis absorbance data
osteo <- read_xlsx("project/data/differentiation/osteogenesis_processed.xlsx")

# DNA concentration data
dna <- read_xlsx("project/data/differentiation/DNA_cell_lysates.xlsx",
                 col_types = c("text", "text", "numeric"))


# Pivot to long format -------------------------------------------

osteo_long <- osteo |>
  pivot_longer(cols = c("0", "8", "8osteo"),
               names_to = "day",
               values_to = "absorption") |>
  mutate(day = as.factor(day))


# Statistical tests ----------------------------------------------

# Clone A vs clone B at day 0, to quantify the effect of clone at
# baseline
osteo_long_day0 <- osteo_long |>
  filter(day == "0")

wilcox_AvB_day0 <- wilcox.test(absorption ~ clone,
                               data = osteo_long_day0,
                               exact = FALSE)

# The Wilcoxon rank-sum test results show no significant difference in
# absorption between clone A and clone B at day 0 (p = 0.077),
# suggesting that the two clones have similar baseline absorption
# levels before the introduction of osteogenic medium. However it
# seems as this is due to low sample size (n = 3). Delta ~= 0.05,
# which is quite large for this data. Thus it seems that the baseline
# absorption is different between the clones, but the statistical
# power to detect this difference is very low.


# Day 8 vs day 0 within each clone
wilcox_A_day8v0 <- wilcox.test(
  osteo_long$absorption[osteo_long$clone == "A" & osteo_long$day == "0"],
  osteo_long$absorption[osteo_long$clone == "A" & osteo_long$day == "8"],
  exact = FALSE)

wilcox_B_day8v0 <- wilcox.test(
  osteo_long$absorption[osteo_long$clone == "B" & osteo_long$day == "0"],
  osteo_long$absorption[osteo_long$clone == "B" & osteo_long$day == "8"],
  exact = FALSE)

# The Wilcoxon rank-sum test results show that the differences in
# absorption between day 0 and day 8 for both clones are not
# statistically significant (clone A: p = 0.077; clone B: p = 0.1).
# However this is also likely due to small sample size, as the data
# does show an increase in absorption (delta(A) ~ 0.17,
# delta(B) ~ 0.2).


# Wilcoxon test on the change in absorption between day 8 and day 0
# between clones
day_change <- osteo_long |>
  group_by(clone, replicate) |>
  summarise(change = absorption[day == "8"] - absorption[day == "0"],
            .groups = "drop")

wilcox_AvB_0to8 <- wilcox.test(change ~ clone, data = day_change)

# The Wilcoxon rank-sum test results show no significant difference in
# the change in absorption between clone A and clone B (p = 0.1),
# suggesting that the overall change in ALP activity from day 0 to
# day 8 is not significantly different between the two clones. delta
# for B looks slightly higher than delta for A, but we will assume
# statistical non-significance for ease of the next step.


# Test change in absorption between day 8 with osteogenic medium
# and day 8 control between clones
osteo_change <- osteo_long |>
  group_by(clone, replicate) |>
  summarise(change = absorption[day == "8osteo"] - absorption[day == "8"],
            .groups = "drop")

wilcox_osteo_change <- wilcox.test(change ~ clone, data = osteo_change)

# The Wilcoxon rank-sum test shows no significant difference between
# the change in the two clones' absorption values between day 8 with
# osteogenic medium and day 8 control (p = 0.1). However, looking at
# the data, the difference between the clones is about 0.196, which
# is quite large.


# Figure for report ----------------------------------------------

# Rescale DNA axis to fit absorbance axis on a dual-axis plot
max_abs <- max(osteo_long$absorption, na.rm = TRUE)
max_dna <- max(dna$conc, na.rm = TRUE)
scale_factor <- max_abs / max_dna

# Relabel clones for nicer facet titles
osteo_long$clone <- factor(osteo_long$clone,
                           levels = c("A", "B"),
                           labels = c("Clone A", "Clone B"))

dna$clone <- factor(dna$clone,
                    levels = c("A", "B"),
                    labels = c("Clone A", "Clone B"))

osteo_plot <- ggplot(osteo_long, aes(x = day)) +
  geom_col(data = dna,
           aes(y = conc * scale_factor, fill = "DNA concentration"),
           position = "dodge", alpha = 0.3, width = 0.6) +
  geom_line(aes(y = absorption, colour = "Absorbance", group = replicate),
            alpha = 0.4, linewidth = 0.7) +
  geom_point(aes(y = absorption, colour = "Absorbance"),
             size = 1) +
  scale_colour_manual(
    name = NULL,
    values = c("Absorbance" = "black"),
    guide = guide_legend(override.aes = list(shape = 16, linetype = 1))) +
  scale_fill_manual(
    name = NULL,
    values = c("DNA concentration" = "grey50")) +
  labs(
    x = "Day",
    y = "Absorbance 405nm (AU)") +
  scale_x_discrete(labels = c(
    "0" = "0",
    "8" = "8",
    "8osteo" = "8 in osteogenic\nmedium")) +
  scale_y_continuous(
    n.breaks = 10,
    sec.axis = sec_axis(~ . / scale_factor,
                        name = expression("[DNA] (µg mL"^-1 * ")"))) +
  facet_wrap2(~ clone, strip = strip_themed(
    background_x = list(
      element_rect(fill = "#f7766f"),
      element_rect(fill = "#0dc1c5")),
    text_x = list(
      element_text(colour = "black", face = "bold"),
      element_text(colour = "black", face = "bold")))) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(colour = "#454644"))
osteo_plot


# Save plot ------------------------------------------------------

ggsave("project/figures/osteo_plot.png",
       plot = osteo_plot,
       width = 6, height = 4, units = "in",
       dpi = 300)


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

# Mangiafico S (2024). _rcompanion: Functions to Support Extension
# Education Program Evaluation_. R package version 2.4.36,
# <https://CRAN.R-project.org/package=rcompanion>.

# Dinno A (2024). _dunn.test: Dunn's Test of Multiple Comparisons
# Using Rank Sums_. R package version 1.3.6,
# <https://CRAN.R-project.org/package=dunn.test>.

# van den Brand T (2024). _ggh4x: Hacks for 'ggplot2'_. R package
# version 0.2.8, <https://CRAN.R-project.org/package=ggh4x>.

# Wickham H, Bryan J (2023). _readxl: Read Excel Files_. R package
# version 1.4.3, <https://CRAN.R-project.org/package=readxl>.