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
tid_mapping <- read_csv("project/data/TID_to_trackingid.csv",
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
  arrange(clone, replicate, lineage.id, tracking.id, frame) |>
  group_by(clone, replicate, lineage.id, tracking.id) |>
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
    mean.speed = mean(instantaneous.velocity, na.rm = TRUE),
    .groups = "drop") # Drop groups after summarising
#
# 3aii. Remove mean.speed NaN objects
livecyte_pretty <- livecyte_pretty |> 
  filter(!is.na(mean.speed))
#
# 3b. Plot the data to determine distributions of common elements: we know that the manual dataset only contains cells, so we'll use that to figure out the parameters to filter the livecyte dataset.
#
# 3bi. Semi-join the livecyte_pretty dataset to a new dataset that contains only the tracking.ids that are present in the manual dataset
livecyte_filtered <- livecyte_pretty %>%
  filter(
    tracking.id %in% tid_mapping$tracking.id,
    (clone == "cloneA" & replicate == 1) |
      (clone == "cloneB" & replicate == 2)
  )
#
# 3bii. Plot the distribution of total_path_length, final_displacement, and mean_speed for the filtered dataset
total_path_length <-ggplot(livecyte_filtered, aes(x = total_path_length)) +
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
# 3c. Discover the thresholds for filtering out debris from the livecyte dataset based on the distributions observed in the manual dataset. We can use the minimum values of track.length, euclidean.distance, and mean.speed from the manual dataset as thresholds for filtering the livecyte dataset.
#
# 3ci. Obtain summary statistics for the livecyte dataset
livecyte_stats <- livecyte_filtered |> 
  group_by(clone, replicate) |> 
  summarise(
    # Use total_path_length instead of track.length
    total_path_length_median = median(total_path_length, na.rm = TRUE),
    total_path_length_q1 = quantile(total_path_length, 0.25, na.rm = TRUE),
    total_path_length_q3 = quantile(total_path_length, 0.75, na.rm = TRUE),
    total_path_length_iqr = total_path_length_q3 - total_path_length_q1,
    total_path_length_min = min(total_path_length, na.rm = TRUE),
    total_path_length_max = max(total_path_length, na.rm = TRUE),
    
    # Use final_displacement instead of euclidean.distance
    final_displacement_median = median(final_displacement, na.rm = TRUE),
    final_displacement_q1 = quantile(final_displacement, 0.25, na.rm = TRUE),
    final_displacement_q3 = quantile(final_displacement, 0.75, na.rm = TRUE),
    final_displacement_min = min(final_displacement, na.rm = TRUE),
    final_displacement_max = max(final_displacement, na.rm = TRUE),
    
    mean.speed_median = median(mean.speed, na.rm = TRUE),
    mean.speed_q1 = quantile(mean.speed, 0.25, na.rm = TRUE),
    mean.speed_q3 = quantile(mean.speed, 0.75, na.rm = TRUE),
    mean.speed_min = min(mean.speed, na.rm = TRUE),
    mean.speed_max = max(mean.speed, na.rm = TRUE),
    .groups = "drop"
  )
#
# 3cii. Obtain the 5th and 95th percentiles for total_path_length, final_displacement, and mean.speed for the whole livecyte dataset
#
livecyte_percentiles <- livecyte_filtered |> 
  group_by(clone, replicate) |> 
  summarise(
    total_path_length_5th = quantile(total_path_length, 0.05, na.rm = TRUE),
    total_path_length_95th = quantile(total_path_length, 0.95, na.rm = TRUE),
    
    final_displacement_5th = quantile(final_displacement, 0.05, na.rm = TRUE),
    final_displacement_95th = quantile(final_displacement, 0.95, na.rm = TRUE),
    
    mean.speed_5th = quantile(mean.speed, 0.05, na.rm = TRUE),
    mean.speed_95th = quantile(mean.speed, 0.95, na.rm = TRUE),
    .groups = "drop"
  )
