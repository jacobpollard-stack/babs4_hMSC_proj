# ============================================================
# Livecyte Data Analysis
# - Data assembling, tidying, and debris filtering
# - DRAFT 2
# ============================================================
#
# 1. Load libraries ---------------------------------------------
#
library(tidyverse)
library(readxl)
library(kSamples)
#
# 2. Load data --------------------------------------------------
#
# 2a. Load the livecyte automatic dataset
#
livecyte <- read_tsv("project/data/movement_morphology/livecyte_data.tsv",
                     col_types = cols(
                       clone = col_factor(),
                       replicate = col_factor(),
                       tracking.id = col_factor(),
                       lineage.id = col_factor()
                     )) # Ensure that categorical variables are read as factors
#
# 2b. Load the manual dataset
#
manual <- read_tsv("project/data/movement_morphology/manual_data.tsv",
                   col_types = cols(
                     track.present.at.start.or.end = col_factor(),
                     track.never.divides = col_factor(),
                     cell.line = col_factor()
                   ))
#
# 2c. Load my TID vs tracking.id dataset; created manually by matching the TIDs from the manual dataset to the tracking.ids from the livecyte dataset, using other parameters.
tid_mapping <- read_csv("project/data/aligning/TID_to_trackingid.csv",
                        col_types = cols(
                          clone = col_factor(),
                          replicate = col_factor(),
                          tracking.id = col_factor(),
                          TID = col_factor(),
                          LID = col_factor(),
                          cell.line = col_factor()
                        ))
#
# 3. Determining what's debris and what's cell ------------------
#
# 3ai. Collapse all the data into one nice table
#
livecyte_pretty <- livecyte |> 
  arrange(clone, replicate, tracking.id, frame) |>
  group_by(clone, replicate, tracking.id) |>
  summarise(
    n_frames = n(),
    start_frame = first(frame),
    total_path_length = last(track.length),
    final_displacement = last(displacement), # As displacement is cumulative
    volume = mean(volume), # Mean over all frames for a single tracking.id within a clone and replicate
    radius = mean(radius),
    sphericity = mean(sphericity),
    length.to.width.ratio = mean(length.to.width),
    dry.mass = mean(dry.mass),
    mean.speed = mean(instantaneous.velocity, na.rm = TRUE)*60, # Convert from µm/sec to µm/min
    .groups = "drop") # Drop groups after summarising
#
# 3aii. Remove mean.speed NaN objects
livecyte_pretty <- livecyte_pretty |> 
  filter(!is.na(mean.speed))
#
# 3b. Plot the data to determine distributions of common elements: we know that the manual dataset only contains cells, so we'll use that to figure out the parameters to filter the livecyte dataset.
#
# 3bi. Filter livecyte_pretty so that it only contains tracking.ids found in tid_mapping, that have the respective clone and replicate.
livecyte_filtered <- livecyte_pretty |> 
  inner_join(tid_mapping, by = c("clone", "replicate", "tracking.id"))
#
# 3bii. Plot the distribution of total_path_length, final_displacement, and mean_speed for the filtered dataset
total_path_length <- ggplot(livecyte_filtered, aes(x = total_path_length)) +
  facet_wrap(~ clone) +
  geom_histogram(binwidth = 75)
final_displacement <- ggplot(livecyte_filtered, aes(x = final_displacement)) +
  facet_wrap(~ clone) +
  geom_histogram(binwidth = 50)
mean.speed_livecyte <- ggplot(livecyte_filtered, aes(x = mean.speed)) +
  facet_wrap(~ clone) +
  geom_histogram(binwidth = 0.001)
# 3biii. Plot the distribution of track.length, euclidean.distance, and mean.speed for the manual dataset
track.length <- ggplot(manual, aes(x = track.length)) +
  facet_wrap(~ cell.line) +
  geom_histogram(binwidth = 150)
euclidean.distance <- ggplot(manual, aes(x = euclidean.distance)) +
  facet_wrap(~ cell.line) +
  geom_histogram(binwidth = 50)
mean.speed_manual <- ggplot(manual, aes(x = mean.speed)) +
  facet_wrap(~ cell.line) +
  geom_histogram(binwidth = 0.05)
# 
# 3biv. Overlay distributions to see if they are similar; if they are, we can use the same thresholds for filtering out debris from the livecyte dataset.
# euclidean.distance and final_displacement
ggplot() +
  geom_histogram(data = manual, aes(x = euclidean.distance), fill = "blue", alpha = 0.5, binwidth = 30) +
  geom_histogram(data = livecyte_filtered, aes(x = final_displacement), fill = "red", alpha = 0.5, binwidth = 30) +
  facet_wrap(~ cell.line)
#
# track.length and total_path_length
ggplot() +
  geom_histogram(data = manual, aes(x = track.length), fill = "blue", alpha = 0.5, binwidth = 30) +
  geom_histogram(data = livecyte_filtered, aes(x = total_path_length), fill = "red", alpha = 0.5, binwidth = 30) +
  facet_wrap(~ cell.line)

# mean.speed
ggplot() +
  geom_histogram(data = manual, aes(x = mean.speed), fill = "blue", alpha = 0.5, binwidth = 0.01) +
  geom_histogram(data = livecyte_filtered, aes(x = mean.speed), fill = "red", alpha = 0.5, binwidth = 0.01) +
  facet_wrap(~ cell.line)

# Clone A appears more likely for the automated Livecyte tracking to lose cells, leading to a dramatically different track length + mean.speed distributions for A1 than B2

# We see clearly that the mean speed and track length distributions are very different between the manual and livecyte datasets, suggesting that the mean speed values in the livecyte dataset are not reliable and should not be used for filtering out debris. The euclidean distance distribution is more similar between the manual and livecyte datasets, suggesting that we can use the same thresholds for filtering out debris from the livecyte dataset based on this parameter. However, we should be cautious when interpreting the results, as there may still be some differences between the two datasets that could affect the analysis.

# 3bv. Perform a KS and  Wilcoxon test to compare the distributions of euclidean distance in manual and final displacement in livecyte_filtered.
ks_test_euc <- ks.test(manual$euclidean.distance, livecyte_filtered$final_displacement)
wilcox_test_euc <- wilcox.test(manual$euclidean.distance, livecyte_filtered$final_displacement)
print(ks_test_euc)
print(wilcox_test_euc)
# KS: D = 0.074627, p-value = 0.6319
# Wilcoxon: W = 20184, p-value = 0.9426
# The KS and Wilcoxon tests suggest that there is no significant difference between the distributions of euclidean distance in the manual dataset and final displacement in the livecyte_filtered dataset, meaning that we can use this parameter for filtering debris from the livecyte dataset.

# 3bvi. Perform a KS and Wilcoxon test to compare the distributions of track.length in manual and total_path_length in livecyte_filtered.
ks_test_track <- ks.test(manual$track.length, livecyte_filtered$total_path_length)
wilcox_test_track <- wilcox.test(manual$track.length, livecyte_filtered$total_path_length)
print(ks_test_track)
print(wilcox_test_track)
# KS: D = 0.30341, p-value = 1.928e-08
# Wilcoxon: W = 22746, p-value = 0.02262
# The KS and Wilcoxon tests suggest that there is a significant difference between the distributions of track.length in the manual dataset and total_path_length in the livecyte_filtered dataset, meaning that we should not use this parameter for filtering debris from the livecyte dataset.
#
# 4. Test whether the distributions of final displacement are roughly homogenous across the different replicates in a single clone.
# We will use an Anderson-Darling k-sample and Kruskal-Wallis tests.
# 4a. Split livecyte_pretty into two datasets, one for each clone.
livecyte_A <- livecyte_pretty |>
  filter(clone == "cloneA")
livecyte_B <- livecyte_pretty |>
  filter(clone == "cloneB")
# 4bi. Perform Anderson-Darling (AD) k-sample test for each new dataset.
# Clone A
ad_test_A <- do.call(ad.test, split(livecyte_A$final_displacement, livecyte_A$replicate))
# Clone B
ad_test_B <- do.call(ad.test, split(livecyte_B$final_displacement, livecyte_B$replicate))
# Print results
print(ad_test_A)
print(ad_test_B)
# The test shows that there is not sufficient evidence to suggest that the distribution of the replicates within clone A are different (p=0.064094), however there is sufficient evidence to suggest that the distribution of the replicates within clone B are different.
# 4bii. Perform pairwise tests on clone B
pairwise.wilcox.test(livecyte_B$final_displacement, 
                     livecyte_B$replicate,
                     p.adjust.method = "BH")
#
#   1       2      
# 2 0.39    -      
# 3 2.1e-11 2.0e-12
# Replicate 3 seems to be significantly different from reps 1 and 2; we will plot this using an empirical cumulative distribution plot.
ggplot(livecyte_B, aes(x = final_displacement, colour = replicate)) +
  stat_ecdf()
# The plot shows that replicate 3 is significantly higher in final displacement than replicates 1 and 2, which are more similar to each other. This suggests that replicate 3 may contain more debris or tracking errors, which could be affecting the distribution of final displacement in clone B.
#
# 4c. Discover whether replicate 3 is different because of a higher proportion of debris, or because the cells in that replicate are simply more motile; ie. is the difference real?
# 4ci. Create a dataset containing the objects in livecyte_pretty that were not matched to a real cell
livecyte_debris <- livecyte_pretty |>
  anti_join(tid_mapping, by = c("clone", "replicate", "tracking.id")) |>
  filter(clone == "cloneB")
# 4cii. Perform statistical tests to compare the distribution of final displacement in the debris dataset to the distribution of final displacement in the matched dataset (livecyte_B).
# AD test
ad_debris_B <- do.call(ad.test, 
                       split(livecyte_debris$final_displacement, 
                             livecyte_debris$replicate))
print(ad_debris_B)
#
#  AD   T.AD  asympt. P-value
# version 1: 28.888 24.982       2.5318e-15
# version 2: 28.900 24.987       2.2875e-15
#
# Wilcoxon test
wilcoxon_debris_B <- pairwise.wilcox.test(livecyte_debris$final_displacement,
                     livecyte_debris$replicate,
                     p.adjust.method = "BH")
print(wilcoxon_debris_B)
#
#   1       2      
# 2 0.079   -      
# 3 4.3e-11 3.8e-10
#
# The tests suggest that the distribution of debris in clone B replicate 3 is significantly different from the distribution of debris in replicates 1 and 2, which are more similar to each other. This suggests that replicate 3 may contain more debris or tracking errors, which could be affecting the distribution of final displacement in clone B.
# Thus, we will filter out replicate 3 from clone B in our livecyte dataset, as it appears to contain a higher proportion of debris or tracking errors, which could be affecting the analysis of cell motility in this clone.
#
# 4ciii.
livecyte_pretty <- livecyte_pretty |>
  filter(!(clone == "cloneB" & replicate == "3"))
#
#
# This is all very complicated and creating more problems than it's solving. Cells and debris overlap significantly across all parameters, making accurate filtering borderline impossible. I'm sure a better statistician could do it, but not I.