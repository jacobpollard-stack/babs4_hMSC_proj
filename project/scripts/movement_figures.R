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
#
# 2. Load data
#
df <- read.delim("project/data/movement_morphology/livecyte_collapsed_filtered.tsv")
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

