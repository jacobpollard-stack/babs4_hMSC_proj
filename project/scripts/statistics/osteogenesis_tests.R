# ==============================================================================
# Livecyte Data Analysis
# - Statistical analysis:
# - Osteogenesis differentiation statistical tests
# =============================================================================
#
# 1. Load libraries
library(tidyverse)
library(rcompanion)
library(dunn.test)
#
# 2. Load data
#
osteo  <- read_xlsx("project/data/differentiation/osteogenesis_processed.xlsx")
#
# 3. Pivot longer
#
osteo_long <- osteo %>%
  pivot_longer(cols = c("0", "8", "8osteo"), names_to = "day", values_to = "absorption") %>%
  mutate(day = as.factor(day))
#
# 4. Statistical tests
#
# 4a. Scheirer–Ray–Hare test on day 8 data
#
# 4ai. Subset data for day 8 and 8osteo to quantify the effect of clone and osteogenic medium
#
osteo_day8 <- osteo_long |>
  filter(day %in% c("8", "8osteo"))
#
# 4aii. Perform test
#
srh_day8 <- scheirerRayHare(absorption ~ clone + day, data = osteo_day8)
#
# 4b. Post-hoc Dunn test
#
# 4bi. Create a combined factor for clone and day
#
osteo_day8 <- osteo_day8 |>
  mutate(group = interaction(clone, day))

dunn_day8 <- dunn.test(osteo_day8$absorption, osteo_day8$group, method = "bh")
#
# Col Mean-│
# Row Mean │        A.8   A.8osteo        B.8
# ─────────┼─────────────────────────────────
# A.8osteo │  -2.038098
#          │     0.0415*
#          │
# B.8      │   1.019049   3.057147
#          │     0.1541     0.0067**
#          │
# B.8osteo │  -1.019049   1.019049  -2.038098
#               0.1849     0.2311     0.0623
#
# The Dunn test results show a significant difference between A.8 and A.8osteo (p = 0.0415) and between A.8osteo and B.8 (p = 0.0067); this is not important. No significant difference exists between the other groups. This suggests that the osteogenic medium has a significant effect on absorption for clone A, but not for clone B (0.0623). In addition, at alpha=0.05, clone A and B did not exhibit significant differences in absorption, regardless of the presence of osteogenic medium. However, the absolute difference in absorption between B.8 and B.8osteo is quite large (~0.282 vs ~0.453), and the p-value is close to the significance threshold, suggesting that there may be a trend towards significance that could be worth further investigation with a larger number of technical replicates.
#
# 5. Clone A vs clone B at day 0
#
# 5a. Subset data for day 0 to quantify the effect of clone at baseline
#
osteo_day0 <- osteo_long |>
  filter(day == "0")
#
# 5b. Perform wilcoxon rank-sum test
#
wilcox_AB_day0 <- wilcox.test(absorption ~ clone, data = osteo_day0)
#
# The Wilcoxon rank-sum test results show no significant difference in absorption between clone A and clone B at day 0 (p = 0.0765225), suggesting that the two clones have similar baseline absorption levels before the introduction of osteogenic medium. 
#
# 6. Day 8 vs day 0 within each clone
#
# 6a. Clone A
#
wilcox_A_day <- wilcox.test(
  osteo_long$absorption[osteo_long$clone == "A" & osteo_long$day == "0"],
  osteo_long$absorption[osteo_long$clone == "A" & osteo_long$day == "8"]
)
#
# 6b. Clone B
#
wilcox_B_day <- wilcox.test(
  osteo_long$absorption[osteo_long$clone == "B" & osteo_long$day == "0"],
  osteo_long$absorption[osteo_long$clone == "B" & osteo_long$day == "8"]
)
#
# The Wilcoxon rank-sum test results show that the differences in absorption between day 0 and day 8 for both clones are not statistically significant (clone A: p = 0.0765225; clone B: p = 0.1). This suggests that ALP  activity did not change significantly over time for either clone in the absence of osteogenic medium, indicating that the observed differences at day 8 may be primarily driven by the presence of osteogenic medium rather than time alone.
#
# 7. Although the Wilcoxon rank-sum test results show no significant differences in absorption between day 0 and day 8 for either clone, it may still be informative to compare the change in absorption over time between the two clones to see if there are any trends that could be worth further investigation with a larger number of technical replicates.
#
# 7a. Wilcoxon test on the change in absorption between day 8 and day between clones
#
# 7ai. Calculate the change in absorption between day 8 and day 0 for each clone
#
day_change <- osteo_long %>%
  group_by(clone, replicate) %>%
  summarise(change = absorption[day == "8"] - absorption[day == "0"])
#
# 7aii. Perform Wilcoxon test on the change in absorption
#
wilcox_day_change <- wilcox.test(change ~ clone, data = day_change)
#
# The Wilcoxon rank-sum test results show no significant difference in the change in absorption between clone A and clone B (p = 0.1), suggesting that the overall change in ALP activity from day 0 to day 8 is not significantly different between the two clones.
#
# 7b. Test change in absorption between day 8 with osteogenic medium and day 8 control between clones
#
# 7bi. Calculate the change in absorption between day 8 with osteogenic medium and day 8 control for each clone
#
osteo_change <- osteo_long %>%
  group_by(clone, replicate) %>%
  summarise(change = absorption[day == "8osteo"] - absorption[day == "8"])
#
# 7bii. Perform Wilcoxon test on the change in absorption between clone A and clone B
#
wilcox_osteo_change <- wilcox.test(change ~ clone, data = osteo_change)
#
# The Wilcoxon rank-sun test shows no significant difference between the change in the two clone' absorption values between day 8 with osteogenic medium and day 8 control (p = 0.1), suggesting that the effect of osteogenic medium on ALP activity is not significantly different between the two clones.
#
# 8. Summary of results
#
# - The Scheirer–Ray–Hare test revealed a significant effect of osteogenic medium on absorption for clone A, but not for clone B.
# - The Dunn post-hoc test showed significant differences between A.8 and A.8osteo, and between A.8osteo and B.8, but not between the other groups. However, the difference between B.8 and B.8osteo was close to significance, suggesting a potential trend that may warrant further investigation.
# - The Wilcoxon rank-sum tests showed no significant differences in absorption between clone A and clone B at day 0, and no significant changes in absorption over time for either clone in the absence of osteogenic medium. Additionally, there were no significant differences in the change in absorption between day 8 and day 0, or between day 8 with osteogenic medium and day 8 control, between the two clones.
#
# Overall, the statistical power of only having 3 technical replicates means that the results should be interpreted with caution, and further investigation with a larger number of replicates is necessary to draw more definitive conclusions about the effects of osteogenic medium and mesenchymal clone on ALP activity.