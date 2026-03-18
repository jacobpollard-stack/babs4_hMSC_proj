# ============================================================
# Livecyte Data Analysis
# - Data assembling, tidying, and debris filtering
# - DRAFT 1
# ============================================================
#
# 1. Load libraries ---------------------------------------------
#
library(tidyverse)
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
# 3. Filter out debris from livecyte dataset and check quality of
# data ----------------------------------------------------------

# 3a. Collapse all the data into one nice table
#
livecyte_pretty <- livecyte |> 
  arrange(clone, replicate, lineage.id, tracking.id, frame) |>
  group_by(clone, replicate, lineage.id, tracking.id) |>
  summarise(
    n_frames = n(),
    start_frame = first(frame),
    total_path_length = last(track.length),
    final_displacement = last(displacement), # As displacement is cumulative
    calculated_displacement = sqrt((last(position.x) - first(position.x))^2 +
                                     (last(position.y) - first(position.y))^2), # Calculate displacement from start to end position to double check that the final displacement is correctly calculated by Livecyte
    volume = mean(volume),
    radius = mean(radius),
    sphericity = mean(sphericity),
    length.to.width.ratio = mean(length.to.width),
    dry.mass = mean(dry.mass),
    mean.speed = mean(instantaneous.velocity),
    .groups = "drop") # Drop groups after summarising
#
# 3b. Check whether the final displacement calculated from the start and end positions matches the cumulative displacement provided by Livecyte. Also check whether lineage.id and tracking.id diverge at any point
#
livecyte_pretty <- livecyte_pretty |>
  mutate(displacement_equal = final_displacement == calculated_displacement)
livecyte_pretty <- livecyte_pretty |>
  mutate(lineage_equal = lineage.id == tracking.id)
# We can see that there are some discrepancies between the final displacement and the calculated displacement, which suggests that there may be some errors in the tracking.
# Therefore we will perform a Wilcoxon test to see if there is a significant difference between the final displacement and the calculated displacement. If there is a significant difference, we will need to investigate further to see if there are any systematic errors in the tracking.
wilcox.test(livecyte_pretty$final_displacement, livecyte_pretty$calculated_displacement, paired = TRUE)
#
# lineage.id and tracking.id are equal for all objects, which suggests that mitosis did not occur during the experiment. However, we can see cell division happening in the video we collected. Thus, the Livecyte data must be wrong or missing some information.
#
# 3c. Filter out debris based on common features of debris, such as:
# - Objects that appear mid-experiment but are not daughter or root cells (left-censored objects), already accounted for
# - Objects with very short path lengths (less than 10 micrometers)
# - Objects with very low displacement (less than 5 micrometers)
# - Objects with short lifespans (less than 5 frames)
# - High aspect ratio objects, like hairs or scratches (length.to.width.ratio >= 20)
# - Volumes outside of the normal range for cells (600-8500 cubic micrometers)
# - Dry mass that deviates highly outside of normal range (150-2100 picograms (Anconelli et al. (2025))
#
cells_livecyte <- livecyte_pretty |>
  filter(
    total_path_length >= 10,
    final_displacement >= 5,
    n_frames >= 5,
    length.to.width.ratio < 20,
    600 < volume & volume < 8500,
    150 < dry.mass & dry.mass < 2100,
    sphericity > 0.12
  )
#
# Create plots to visualise the distribution of the features before and after filtering to check that the filtering process is working as expected
#
drymass_cloneA <- ggplot(livecyte_pretty, aes(x = dry.mass)) +
  geom_histogram(binwidth = 50, fill = "lightblue", color = "black") +
  facet_wrap(~ clone) +
  labs(title = "Distribution of Dry Mass Before Filtering", x = "Dry Mass (pg)", y = "Count")
drymass_cloneA
drymass_cloneA_filtered <- ggplot(cells_livecyte, aes(x = dry.mass)) +
  geom_histogram(binwidth = 50, fill = "lightblue", color = "black") +
  facet_wrap(~ clone) +
  labs(title = "Distribution of Dry Mass After Filtering", x = "Dry Mass (pg)", y = "Count")
drymass_cloneA_filtered
sphericity_cloneA <- ggplot(livecyte_pretty, aes(x = sphericity)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black") +
  facet_wrap(~ clone) +
  labs(title = "Distribution of Sphericity Before Filtering", x = "Sphericity", y = "Count")
sphericity_cloneA
sphericity_cloneA_filtered <- ggplot(cells_livecyte, aes(x = sphericity)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black") +
  facet_wrap(~ clone) +
  labs(title = "Distribution of Sphericity After Filtering", x = "Sphericity", y = "Count")
sphericity_cloneA_filtered
length_to_width_cloneA <- ggplot(livecyte_pretty, aes(x = length.to.width.ratio)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  facet_wrap(~ clone) +
  labs(title = "Distribution of Length to Width Ratio Before Filtering", x = "Length to Width Ratio", y = "Count")
length_to_width_cloneA
length_to_width_cloneA_filtered <- ggplot(cells_livecyte, aes(x = length.to.width.ratio)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  facet_wrap(~ clone) +
  labs(title = "Distribution of Length to Width Ratio After Filtering", x = "Length to Width Ratio", y = "Count")
length_to_width_cloneA_filtered
#
manual_speed_cloneA <- ggplot(manual, aes(x = mean.speed)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
  facet_wrap(~ cell.line) +
  labs(title = "Distribution of Instantaneous Velocity in Manual Dataset", x = "Instantaneous Velocity (µm/min)", y = "Count")
manual_speed_cloneA
filtered_livecyte_speed <- ggplot(cells_livecyte, aes(x = mean.speed)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
  facet_wrap(~ clone) +
  labs(title = "Distribution of Mean Speed in Filtered Livecyte Dataset", x = "Mean Speed (µm/min)", y = "Count")
filtered_livecyte_speed
#
#
#
# Okay, the debris filtering is defo a good idea, but it needs to be faceted by replicate/clone!