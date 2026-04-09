# =============================================================================
# Livecyte Data Analysis
# - Data assembling, tidying, and debris filtering
# =============================================================================
#
# 1. Load libraries -----------------------------------------------------------
#
library(tidyverse)
library(readxl)
#
# 2. Load data ----------------------------------------------------------------
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
tid_mapping <- read_csv("project/data/movement_morphology/TID_to_trackingid.csv",
                        col_types = cols(
                          clone = col_factor(),
                          replicate = col_factor(),
                          tracking.id = col_factor(),
                          TID = col_factor(),
                          LID = col_factor(),
                          cell.line = col_factor()
                        ))
#
# 3. Determining what's debris and what's cell: -------------------------------
# The livecyte data contains a lot of debris, which overlaps heavily with cells. Thus, we will need to use the  manual data, aligned with cloneA rep1 and cloneB rep2 to determine which tracks in the livecyte data are debris and which are cells. This will further allow us to look at how the debris is distributed, and whether it is more likely to be tracked in certain clones or replicates.
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
    mean.speed = mean(instantaneous.velocity, na.rm = TRUE) *  60, # Convert from microm/sec to microm/h
    .groups = "drop") # Drop groups after summarising
#
# 3aii. Remove mean.speed NaN objects, which are likely debris, as they have no speed (eg. tracked only for one frame)
livecyte_pretty <- livecyte_pretty |> 
  filter(!is.na(mean.speed))
#
# 3b. Now we have a dataset with one row for tracking.id, we can filter it based on the manual dataset.
# 3bi. Filter livecyte_pretty so that it only contains tracking.ids found in tid_mapping, that have the respective clone and replicate.
livecyte_matched <- livecyte_pretty |> 
  inner_join(tid_mapping, by = c("clone", "replicate", "tracking.id"))
#
# 3bii. Plot the distribution of total_path_length, final_displacement, and mean_speed for the filtered dataset
total_path_length <- ggplot(livecyte_matched, aes(x = total_path_length)) +
  facet_wrap(~ clone) +
  geom_histogram(binwidth = 75)
final_displacement <- ggplot(livecyte_matched, aes(x = final_displacement)) +
  facet_wrap(~ clone) +
  geom_histogram(binwidth = 50)
mean.speed_livecyte <- ggplot(livecyte_matched, aes(x = mean.speed)) +
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
# 3biv. Overlay distributions to see if they are similar
# final_displacement and euclidean.distance
ggplot() +
  geom_histogram(data = manual, aes(x = euclidean.distance), fill = "blue", alpha = 0.5, binwidth = 30) +
  geom_histogram(data = livecyte_matched, aes(x = final_displacement), fill = "red", alpha = 0.5, binwidth = 30) +
  facet_wrap(~ cell.line)
#
# track.length and total_path_length
ggplot() +
  geom_histogram(data = manual, aes(x = track.length), fill = "blue", alpha = 0.5, binwidth = 30) +
  geom_histogram(data = livecyte_matched, aes(x = total_path_length), fill = "red", alpha = 0.5, binwidth = 30) +
  facet_wrap(~ cell.line)
#
# mean.speed
ggplot() +
  geom_histogram(data = manual, aes(x = mean.speed), fill = "blue", alpha = 0.5, binwidth = 0.03) +
  geom_histogram(data = livecyte_matched, aes(x = mean.speed), fill = "red", alpha = 0.5, binwidth = 0.03) +
  facet_wrap(~ cell.line)
#
# Based on the distributions, the metrics for the cells in cloneA seems to be more likely to deviate from the manual dataset. This suggests that the tracking is less accurate for cloneA.
#
# 4. Downstream data filtering ------------------------------------------------
# 4a. Create a list of parameters to filter on
filter_params <- list(
  n_frames = 5,
  total_path_length = 100,
  final_displacement = 50,
  volume = 350,
  radius = 15,
  sphericity = 0.1,
  length.to.width.ratio = 1.3,
  dry.mass = 90,
  mean.speed = 0.1
)
#
# 4b. Filter the dataset based on these parameters
livecyte_final <- livecyte_pretty |> 
  filter(
    n_frames >= filter_params$n_frames &
    total_path_length >= filter_params$total_path_length &
    final_displacement >= filter_params$final_displacement  &
    volume >= filter_params$volume &
    radius >= filter_params$radius &
    sphericity >= filter_params$sphericity &
    length.to.width.ratio >= filter_params$length.to.width.ratio  &
    dry.mass >= filter_params$dry.mass &
    mean.speed >= filter_params$mean.speed
  )
#
# 4c. Change cloneA and cloneB to simply A and B to make figure legends nicer
#
livecyte_final <- livecyte_final |> 
  mutate(clone = recode(clone, "cloneA" = "A", "cloneB" = "B"))
#
# 5. Only keep tracking.ids in the livecyte dataset that are found in our new filtered dataset --------------------------------------------------------------
livecyte <- livecyte_final |>
  select(clone, replicate, tracking.id) |>
  inner_join(livecyte, by = c("clone", "replicate", "tracking.id"))
#
# 6. Save the final filtered dataset ------------------------------------------
write_tsv(livecyte_final, "project/data/movement_morphology/livecyte_collapsed_filtered.tsv")
#
# 7. Save the final filtered dataset with all frames --------------------------
write_tsv(livecyte, "project/data/movement_morphology/livecyte_filtered.tsv")
