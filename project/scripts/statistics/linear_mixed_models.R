# ==============================================================================
# Livecyte Data Analysis
# - Statistical analysis:
# - Linear Mixed Models
# ==============================================================================
#
# We are performing a LMM to account for the statistical non-indepencence of the data for replicates within a clone. 
#
# 1. Load libraries ------------------------------------------------------------
#
library(lme4)
library(tidyverse)
library(lmerTest)
library(MuMIn)
#
# 2. Load data -----------------------------------------------------------------
#
collapsed <- read_tsv('project/data/movement_morphology/livecyte_collapsed_filtered.tsv',
                      col_types = cols(
                        clone = col_factor(),
                        replicate = col_factor(),
                        tracking.id = col_factor(),
                        lineage.id = col_factor()
                      ))
#
livecyte <- read_tsv('project/data/movement_morphology/livecyte_filtered.tsv',
                     col_types = cols(
                       clone = col_factor(),
                       replicate = col_factor(),
                       tracking.id = col_factor(),
                       lineage.id = col_factor()
                     ))
#
# 3. Kruskal-Wallace for independence ------------------------------------------
#
# We are interested in whether each parameter is independent between replicates within a clone.
#
# 3a. What features are we interested in?
#
features <- c("dry.mass", "volume", "radius", "sphericity",
              "length.to.width.ratio", "mean.speed",
              "total_path_length", "final_displacement")
#
# 3b. Create a dataset composed of KW test results incl. effect size
#
results <- expand_grid(clone = unique(collapsed$clone), feature = features) |>
  rowwise() |>
  mutate(
    test = list({
      d <- collapsed[collapsed$clone == clone, ]
      kruskal.test(reformulate("replicate", feature), data = d)
    }),
    H        = test$statistic,
    p.value  = test$p.value,
    n        = sum(collapsed$clone == clone),
    epsilon2 = H / (n - 1)
  ) |>
  ungroup() |>
  mutate(p.adjusted = p.adjust(p.value, method = "BH")) |>
  select(clone, feature, H, p.value, p.adjusted, epsilon2)
#
# clone  feature                    H   p.value p.adjusted epsilon2
# cloneA dry.mass               3.26  0.196       0.262     0.00491 
# cloneA volume                 3.26  0.196       0.262     0.00491 
# cloneA radius                 2.10  0.350       0.431     0.00316 
# cloneA sphericity             0.709 0.701       0.701     0.00107 
# cloneA length.to.width.ratio 13.5   0.00118     0.00630   0.0203  ***
# cloneA mean.speed            22.6   0.0000122   0.000195  0.0341  ***
# cloneA total_path_length      5.19  0.0748      0.121     0.00782 
# cloneA final_displacement     5.17  0.0755      0.121     0.00779 
# cloneB dry.mass               9.19  0.0101      0.0323    0.0150  ***
# cloneB volume                 9.19  0.0101      0.0323    0.0150  ***
# cloneB radius                 7.58  0.0225      0.0515    0.0124  ***
# cloneB sphericity             1.17  0.558       0.595     0.00191 
# cloneB length.to.width.ratio  1.74  0.418       0.478     0.00285 
# cloneB mean.speed            17.6   0.000148    0.00118   0.0289  ***
# cloneB total_path_length      6.88  0.0320      0.0641    0.0113  ***
# cloneB final_displacement     7.66  0.0217      0.0515    0.0125  ***
#
# Effect sizes are generally low, signifying that intereplicate variability only accounts for a small proportion of the variance in the data. 
#
# 4. Linear Mixed Models ------------------------------------------------------
#
# 4a. Fit LMM with clone as a fixed effect and replicate as a random effect, creating a table with all combinations of features and clones.
#
lmm_results <- map_dfr(features, function(f) {
  model <- lmer(reformulate(c("clone", "(1 | clone:replicate)"), response = f),
                data = collapsed)
  s <- summary(model)
  coef_name <- rownames(s$coefficients)[2]
  ci <- confint(model, parm = coef_name, method = "Wald")
  r2 <- r.squaredGLMM(model)
  vc <- as.data.frame(VarCorr(model))
  tibble(
    feature           = f,
    estimate          = s$coefficients[coef_name, "Estimate"],
    se                = s$coefficients[coef_name, "Std. Error"],
    ci_lower          = ci[1],
    ci_upper          = ci[2],
    t_value           = s$coefficients[coef_name, "t value"],
    df                = s$coefficients[coef_name, "df"],
    p_value           = s$coefficients[coef_name, "Pr(>|t|)"],
    marginal_r2       = r2[1, "R2m"],
    conditional_r2    = r2[1, "R2c"],
    replicate_var     = vc$vcov[vc$grp == "clone:replicate"],
    residual_var      = vc$vcov[vc$grp == "Residual"],
    icc               = vc$vcov[vc$grp == "clone:replicate"] / 
      sum(vc$vcov)
  )
})
#
# 4b. Adjust p-values and add stars for significance, then form a human-readable table.
#
lmm_results <- lmm_results |>
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ ""))
#
# 4c. Plot results with significance stars
#
ggplot(lmm_results, aes(x = feature, y = estimate)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) + 
  geom_text(aes(label = sig), vjust = -1) + 
  theme_bw() + 
  labs(title = "Linear Mixed Model Results",
       x = "Feature",
       y = "Estimated Effect of Clone") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
