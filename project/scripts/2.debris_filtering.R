# Data assembling, tidying, and debris filtering -----------------

# ================================================================

# Experimental overview ------------------------------------------

# The automatically tracked livecyte dataset contains a lot of debris,
# which overlaps heavily with cells. We use the manual dataset to look
# at how the cells are actually distributed, then apply filtering
# thresholds to remove debris from the livecyte data.


# Description of data --------------------------------------------

# Three data files are used in this script:
#    livecyte_data.tsv                 — per-frame automatic tracking data
#    livecyte_collapsed_unfiltered.tsv — one row per tracking.id, created
#                                        by 1.tid_mapping.R
#    manual_data.tsv                   — manually tracked cells
#    TID_to_trackingid.csv             — mapping between manual TIDs and
#                                        livecyte tracking.ids, created
#                                        by 1.tid_mapping.R
#
# All data is stored in project/data/movement_morphology/.


# Analysis overview ----------------------------------------------

# This script loads the collapsed livecyte dataset and the manual
# dataset, then overlays their distributions to check whether the
# matched cells have similar metrics. We then apply a set of
# filtering thresholds to remove debris from the livecyte data, and
# save both the collapsed filtered dataset and the per-frame filtered
# dataset for downstream analysis.

# ================================================================

# Packages required ----------------------------------------------

# for data manipulation and visualisation
library(tidyverse)

# for reading Excel files
library(readxl)


# Data import ----------------------------------------------------

# raw livecyte automatic dataset (per-frame)
livecyte <- read_tsv("project/data/movement_morphology/livecyte_data.tsv",
                     col_types = cols(
                       clone = col_factor(),
                       replicate = col_factor(),
                       tracking.id = col_factor(),
                       frame = col_integer()
                     ))

# collapsed livecyte automatic dataset (one row per tracking.id)
livecyte_collapsed_unfiltered <- read_tsv(
  "project/data/movement_morphology/livecyte_collapsed_unfiltered.tsv",
  col_types = cols(
    clone = col_factor(),
    replicate = col_factor(),
    tracking.id = col_factor(),
    lineage.id = col_factor()
  ))

# manual dataset
manual <- read_tsv("project/data/movement_morphology/manual_data.tsv",
                   col_types = cols(
                     track.present.at.start.or.end = col_factor(),
                     track.never.divides = col_factor(),
                     cell.line = col_factor()
                   ))


# Comparing manual and livecyte distributions --------------------

# We use the manual data, aligned with clone A rep 1 and clone B
# rep 2, to determine which tracks in the livecyte data are debris
# and which are cells. This will further allow us to look at how the
# debris is distributed, and whether it is more likely to be tracked
# in certain clones or replicates.

# Load TID to tracking.id mapping, created by matching the TIDs
# from the manual dataset to the tracking.ids from the livecyte
# dataset using euclidean distance plus/minus error bounds
tid_mapping <- read_csv(
  "project/data/movement_morphology/TID_to_trackingid.csv",
  col_types = cols(
    clone = col_factor(),
    replicate = col_factor(),
    tracking.id = col_factor(),
    TID = col_factor(),
    LID = col_factor(),
    cell.line = col_factor()
  ))

# Filter livecyte_collapsed_unfiltered so that it only contains
# tracking.ids found in tid_mapping with the respective clone and
# replicate
livecyte_collapsed_unfiltered_matched <- livecyte_collapsed_unfiltered |>
  inner_join(tid_mapping, by = c("clone", "replicate", "tracking.id"))

# Overlay distributions to check whether the matched livecyte cells
# have similar metric distributions to the manual cells

# final_displacement vs euclidean.distance
ggplot() +
  geom_histogram(data = manual,
                 aes(x = euclidean.distance),
                 fill = "blue", alpha = 0.5, binwidth = 30) +
  geom_histogram(data = livecyte_collapsed_unfiltered_matched,
                 aes(x = final_displacement),
                 fill = "red", alpha = 0.5, binwidth = 30) +
  facet_wrap(~cell.line)

# total_path_length vs track.length
ggplot() +
  geom_histogram(data = manual,
                 aes(x = track.length),
                 fill = "blue", alpha = 0.5, binwidth = 30) +
  geom_histogram(data = livecyte_collapsed_unfiltered_matched,
                 aes(x = total_path_length),
                 fill = "red", alpha = 0.5, binwidth = 30) +
  facet_wrap(~ cell.line)

# mean.speed
ggplot() +
  geom_histogram(data = manual,
                 aes(x = mean.speed),
                 fill = "blue", alpha = 0.5, binwidth = 0.03) +
  geom_histogram(data = livecyte_collapsed_unfiltered_matched,
                 aes(x = mean.speed),
                 fill = "red", alpha = 0.5, binwidth = 0.03) +
  facet_wrap(~ cell.line)

# Based on the distributions, the metrics for the cells in clone A
# seem to be more likely to deviate from the manual dataset. This
# suggests that the tracking is less accurate for clone A. In theory,
# the distributions should be near identical, but the livecyte dataset
# contains a lot of debris and the manual dataset is very small, so we
# won't be able to use it to filter the livecyte data directly.


# Downstream data filtering --------------------------------------

# These are arbitrary parameters chosen based on the distributions of
# the livecyte dataset to remove obvious debris
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

# Apply filters
livecyte_collapsed_filtered <- livecyte_collapsed_unfiltered |>
  filter(
    n_frames >= filter_params$n_frames &
      total_path_length >= filter_params$total_path_length &
      final_displacement >= filter_params$final_displacement &
      volume >= filter_params$volume &
      radius >= filter_params$radius &
      sphericity >= filter_params$sphericity &
      length.to.width.ratio >= filter_params$length.to.width.ratio &
      dry.mass >= filter_params$dry.mass &
      mean.speed >= filter_params$mean.speed
  )


# Filter per-frame data to match --------------------------------

# Only keep tracking.ids in the raw livecyte dataset that are found
# in our new filtered dataset
livecyte_data_filtered <- livecyte |>
  semi_join(livecyte_collapsed_filtered,
            by = c("clone", "replicate", "tracking.id"))


# Save filtered datasets -----------------------------------------

# collapsed (one row per tracking.id)
write_tsv(livecyte_collapsed_filtered,
          "project/data/movement_morphology/livecyte_collapsed_filtered.tsv")

# per-frame (all frames for filtered tracking.ids)
write_tsv(livecyte_data_filtered,
          "project/data/movement_morphology/livecyte_data_filtered.tsv")


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

# Wickham H, Bryan J (2023). _readxl: Read Excel Files_. R package
# version 1.4.3, <https://CRAN.R-project.org/package=readxl>.