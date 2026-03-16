# ============================================================
# Livecyte Data Analysis
# - Data assembling, tidying, and debris filtering
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
# 3. Filter out debris from livecyte dataset --------------------

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
    .groups = "drop") |> # Drop groups after summarising
  mutate(
    lineage_rootcell = start_frame == min(start_frame),
    left_censored = ! lineage_rootcell & start_frame > 1) # Flags objects that appear mid-experiment but are not daughter or root cells ie. debris
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
livecyte_pretty |>
  filter(!lineage_equal)
# lineage.id and tracking.id are equal for all objects, which suggests that mitosis did not occur during the experiment. However, we can see cell division happening in the video we collected. Thus, the Livecyte data must be wrong or missing some information.
#
# 3c. Collapsing the manual dataset
#



# 3d. Initially, we can eliminate debris from colony A, replicate 1 and colony B, replicate 2 by using the manual dataset that only includes cells and no debris.
# 


