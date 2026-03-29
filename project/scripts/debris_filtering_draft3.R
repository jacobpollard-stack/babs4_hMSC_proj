# ============================================================
# Livecyte Data Analysis
# - Data assembling, tidying, and debris filtering
# - DRAFT 3
# ============================================================
#
# 1. Load libraries ------------------------------------------
#
library(tidyverse)
library(readxl)
#
# 2. Load data -----------------------------------------------
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
# 2c. Load my TID vs tracking.id dataset; created manually by matching the TIDs from the manual dataset to the tracking.ids from the livecyte dataset, using euclidean distance plus/minus error bounds.
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
# 3. Determining what's debris and what's cell:
# The livecyte data contains a lot of debris, which overlaps heavily with cells. Thus, we will need to use the  manual data, aligned with cloneA rep1 and cloneB rep2 to determine which tracks in the livecyte data are debris and which are cells.
#
# 3ai. The livecyte dataset contains data for all 'frames' taken over the 4-day period, so we will collapse the frames into a single row per tracking.id.
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
# 3aii. Remove mean.speed NaN objects, which are likely debris, as they have no speed (eg. tracked only for one frame)
livecyte_pretty <- livecyte_pretty |> 
  filter(!is.na(mean.speed))
#
# 3b. Plot the data to determine distributions of tracking perameters. This may allow us to determine which tracks are debris and which are cells.
#
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
#
# mean.speed
ggplot() +
  geom_histogram(data = manual, aes(x = mean.speed), fill = "blue", alpha = 0.5, binwidth = 0.01) +
  geom_histogram(data = livecyte_filtered, aes(x = mean.speed), fill = "red", alpha = 0.5, binwidth = 0.01) +
  facet_wrap(~ cell.line)
#
# Clone A appears more likely for the automated Livecyte tracking to lose cells, leading to a different distributions for A1 than B2.
#
# 3c. Perform Wilcoxon rank sum and Kolmogorov-Smirnov tests to determine if the distributions of the manual and livecyte_filtered datasets are significantly different for each parameter. This will allow us to quantify the accuracy of the livecyte tracking data.
#
# 3ci. euclidean.distance and final_displacement
ks_test_euc <- ks.test(manual$euclidean.distance, livecyte_filtered$final_displacement)
wilcox_test_euc <- wilcox.test(manual$euclidean.distance, livecyte_filtered$final_displacement)
print(ks_test_euc)
print(wilcox_test_euc)
# KS: D = 0.074627, p-value = 0.6319
# Wilcoxon: W = 20184, p-value = 0.9426
#
# 3cii. track.length and total_path_length
ks_test_track <- ks.test(manual$track.length, livecyte_filtered$total_path_length)
wilcox_test_track <- wilcox.test(manual$track.length, livecyte_filtered$total_path_length)
print(ks_test_track)
print(wilcox_test_track)
# KS: D = 0.30341, p-value = 1.928e-08
# Wilcoxon: W = 22746, p-value = 0.02262
#
# 3ciii. mean.speed
ks_test_speed <- ks.test(manual$mean.speed, livecyte_filtered$mean.speed)
wilcox_test_speed <- wilcox.test(manual$mean.speed, livecyte_filtered$mean.speed)
print(ks_test_speed)
print(wilcox_test_speed)
# KS: D = 0.19, p-value = 0.001464
# W = 21994, p-value = 0.08464
#
# Based on the statistical tests, the distributions 