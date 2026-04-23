# Linear mixed models for inter-replicate variability ------------

# ================================================================

# Experimental overview ------------------------------------------

# We are comparing morphological and movement features between two
# clonal cell lines (A and B), each with three biological replicates.
# Because replicates within a clone are not statistically independent,
# we need to account for this using linear mixed models (LMMs).


# Description of data --------------------------------------------

# The filtered livecyte dataset (livecyte_collapsed_filtered.tsv)
# contains one row per tracking.id with columns: clone, replicate,
# tracking.id, and nine numeric features (dry.mass, volume, radius,
# sphericity, length.to.width.ratio, mean.speed, total_path_length,
# final_displacement, mean.thickness). The data is stored in
# project/data/movement_morphology/.


# Analysis overview ----------------------------------------------

# This script first performs Kruskal-Wallis tests to quantify the
# effect of replicate within each clone, and calculates epsilon-squared
# effect sizes to understand how much of the variance is due to
# inter-replicate variability. We then fit LMMs with clone as a fixed
# effect and replicate nested within clone as a random effect, extract
# effect sizes (in SD units), confidence intervals, and ICC values,
# and visualise the results in a forest plot. Results are saved in
# project/figures/.

# ================================================================

# Packages required ----------------------------------------------

# for linear mixed models
library(lme4)

# for data manipulation and plotting
library(tidyverse)

# for p-values in LMMs
library(lmerTest)

# for R-squared in LMMs
library(MuMIn)


# Data import ----------------------------------------------------

livecyte_collapsed_filtered <- read_tsv(
  "project/data/movement_morphology/livecyte_collapsed_filtered.tsv",
  col_types = cols(
    clone = col_factor(),
    replicate = col_factor(),
    tracking.id = col_factor()
  ))


# Define features and standardise --------------------------------

features <- c("dry.mass", "volume", "radius", "sphericity",
              "length.to.width.ratio", "mean.speed",
              "total_path_length", "final_displacement", "mean.thickness")

# Standardise all feature columns so values are in standard deviation
# units, similar to Cohen's d. This allows us to compare effect sizes
# across features on a common scale
livecyte_collapsed_filtered_scaled <- livecyte_collapsed_filtered |>
  mutate(across(all_of(features), ~ scale(.x)[, 1]))


# Kruskal-Wallis testing -----------------------------------------

# We are interested in whether each parameter is independent between
# replicates within a clone. This tells us whether inter-replicate
# variability is a concern that needs to be accounted for in the
# statistical model

results <- expand_grid(
  clone = unique(livecyte_collapsed_filtered_scaled$clone),
  feature = features
) |>
  rowwise() |>
  mutate(
    test = list({
      d <- livecyte_collapsed_filtered[livecyte_collapsed_filtered$clone == clone, ]
      kruskal.test(reformulate("replicate", feature), data = d)
    }),
    H        = test$statistic,
    p.value  = test$p.value,
    n        = sum(livecyte_collapsed_filtered$clone == clone),
    epsilon2 = H / (n - 1)
  ) |>
  ungroup() |>
  mutate(p.adjusted = p.adjust(p.value, method = "BH")) |>
  select(clone, feature, H, p.value, p.adjusted, epsilon2)

# clone  feature                    H   p.value p.adjusted epsilon2
# cloneA dry.mass               3.26  0.196       0.262     0.00491
# cloneA volume                 3.26  0.196       0.262     0.00491
# cloneA radius                 2.10  0.350       0.431     0.00316
# cloneA sphericity             0.709 0.701       0.701     0.00107
# cloneA length.to.width.ratio 13.5   0.00118     0.00630   0.0203
# cloneA mean.speed            22.6   0.0000122   0.000195  0.0341
# cloneA total_path_length      5.19  0.0748      0.121     0.00782
# cloneA final_displacement     5.17  0.0755      0.121     0.00779
# cloneA mean.thickness         2.06  0.356       0.458     0.00311
# cloneB dry.mass               9.19  0.0101      0.0323    0.0150
# cloneB volume                 9.19  0.0101      0.0323    0.0150
# cloneB radius                 7.58  0.0225      0.0515    0.0124
# cloneB sphericity             1.17  0.558       0.595     0.00191
# cloneB length.to.width.ratio  1.74  0.418       0.478     0.00285
# cloneB mean.speed            17.6   0.000148    0.00118   0.0289
# cloneB total_path_length      6.88  0.0320      0.0641    0.0113
# cloneB final_displacement     7.66  0.0217      0.0515    0.0125
# cloneB mean.thickness         0.893 0.640       0.677     0.00146

# Effect sizes are generally low (all epsilon-squared < 0.05),
# signifying that inter-replicate variability only accounts for a
# small proportion of the variance in the data. This supports using
# a random intercept for replicate rather than treating replicates
# as entirely separate groups.


# Linear mixed models --------------------------------------------

# Fit LMM with clone as a fixed effect and replicate nested within
# clone as a random effect. We extract the estimate, standard error,
# confidence interval, t-value, degrees of freedom, p-value,
# marginal and conditional R-squared, replicate variance, residual
# variance, and ICC for each feature
lmm_results <- map_dfr(features, function(f) {
  model <- lmer(
    reformulate(c("clone", "(1 | clone:replicate)"), response = f),
    data = livecyte_collapsed_filtered_scaled
  )
  s <- summary(model)
  coef_name <- rownames(s$coefficients)[2]
  ci <- confint(model, parm = coef_name, method = "Wald")
  r2 <- r.squaredGLMM(model)
  vc <- as.data.frame(VarCorr(model))
  
  tibble(
    feature        = f,
    estimate       = s$coefficients[coef_name, "Estimate"],
    se             = s$coefficients[coef_name, "Std. Error"],
    ci_lower       = ci[1],
    ci_upper       = ci[2],
    t_value        = s$coefficients[coef_name, "t value"],
    df             = s$coefficients[coef_name, "df"],
    p_value        = s$coefficients[coef_name, "Pr(>|t|)"],
    marginal_r2    = r2[1, "R2m"],
    conditional_r2 = r2[1, "R2c"],
    replicate_var  = vc$vcov[vc$grp == "clone:replicate"],
    residual_var   = vc$vcov[vc$grp == "Residual"],
    icc            = vc$vcov[vc$grp == "clone:replicate"] / sum(vc$vcov)
  )
})

# Adjust p-values for multiple comparisons and add significance stars
lmm_results <- lmm_results |>
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    sig = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01  ~ "**",
      p_adjusted < 0.05  ~ "*",
      TRUE               ~ ""
    )
  )


# Figure for report ----------------------------------------------

lmm_results_plot <- ggplot(lmm_results, aes(x = feature, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_text(aes(label = sig), vjust = -1.5) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "black",
             linewidth = 0.1) +
  theme_bw() +
  labs(x = "Feature",
       y = "Effect of Clone (SD units)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(-2, 2))
lmm_results_plot


# Save figure ----------------------------------------------------

units <- "in"
fig_w <- 5
fig_h <- 4
dpi <- 300
device <- "png"

ggsave("project/figures/lmm_results_plot.png",
       plot = lmm_results_plot,
       device = device,
       width = fig_w,
       height = fig_h,
       units = units,
       dpi = dpi)


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

# Bates D, Maechler M, Bolker B, Walker S (2015). "Fitting Linear
# Mixed-Effects Models Using lme4." _Journal of Statistical Software_,
# *67*(1), 1-48. doi:10.18637/jss.v067.i01
# <https://doi.org/10.18637/jss.v067.i01>.

# Kuznetsova A, Brockhoff PB, Christensen RHB (2017). "lmerTest
# Package: Tests in Linear Mixed Effects Models." _Journal of
# Statistical Software_, *82*(13), 1-26. doi:10.18637/jss.v082.i13
# <https://doi.org/10.18637/jss.v082.i13>.

# Barton K (2024). _MuMIn: Multi-Model Inference_. R package version
# 1.48.4, <https://CRAN.R-project.org/package=MuMIn>.