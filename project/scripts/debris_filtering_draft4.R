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
# 3biv. Overlay distributions to see if they are similar
# final_displacement and euclidean.distance
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
  geom_histogram(data = manual, aes(x = mean.speed), fill = "blue", alpha = 0.5, binwidth = 0.03) +
  geom_histogram(data = livecyte_filtered, aes(x = mean.speed), fill = "red", alpha = 0.5, binwidth = 0.03) +
  facet_wrap(~ cell.line)
#
# Based on the distributions, livecyte seems much more likely to track debris for A1 than B2. final_displacement and euclidean.distance seem to be the most reliable metric.
#
# 4. Create a 'debris' dataset and split it into A and B
# 4a.
debris <- livecyte_pretty |> 
  anti_join(tid_mapping, by = c("clone", "replicate", "tracking.id"))
#
# 4b.
debris_A <- debris |> 
  filter(clone == "cloneA" & replicate == "1")
debris_B <- debris |> 
  filter(clone == "cloneB" & replicate == "2")
#
# 4c.
debris_A_other <- debris |> 
  filter(clone == "cloneA" & replicate == "2" | clone == "cloneA" & replicate == "3")
debris_B_other <- debris |> 
  filter(clone == "cloneB" & replicate == "1" | clone == "cloneB" & replicate == "3")
#
# 5. Use Wilcox tests to see whether the debris distributions differ significantly from the debris in A1 and B2.
# A
wilcox.test(debris_A$final_displacement, debris_A_other |> filter(replicate == "2") |> pull(final_displacement))
wilcox.test(debris_A$final_displacement, debris_A_other |> filter(replicate == "3") |> pull(final_displacement))
# B
wilcox.test(debris_B$final_displacement, debris_B_other |> filter(replicate == "1") |> pull(final_displacement))
wilcox.test(debris_B$final_displacement, debris_B_other |> filter(replicate == "3") |> pull(final_displacement))
#
# The debris distribution for B3 is significantly different from the debris in B2 (W = 1948729, p-value = 0.0000000002552), while distributions for the other replicates are not significantly different. As the other debris distributions are relatively homogenous, they should all skew the data roughly in the same way, so we can keep all the debris in the dataset, but we should exclude B3 as it will skew the data more than the other replicates.
#
# 6. Filter the livecyte dataset to exclude B3
livecyte_excl <- livecyte_pretty |> 
  filter(!(clone == "cloneB" & replicate == "3"))
#
# 7. Downstream data filtering
# 7a. Create a list of parameters to filter on
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
# 7b. Filter the dataset based on these parameters
livecyte_final <- livecyte_excl |> 
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
# 8. Save the final filtered dataset
write_tsv(livecyte_final, "project/data/movement_morphology/livecyte_filtered.tsv")
