# Mapping manual tracking IDs to the livecyte dataset ------------

# ================================================================

# Experimental overview ------------------------------------------

# In the experiment, Livecyte used its automatic cell tracking software
# to assign a variety of metrics to each cell track, with a frame taken
# every 23 minutes for 4 days. This tracking is inaccurate and very
# noisy, and the tracking IDs (tracking.id) do not directly
# correspond to the manual tracking IDs (TID) in the manual ImageJ
# dataset. In order to perform a direct comparison between the two
# datasets, we need to match each TID to the corresponding tracking.id.


# Description of data --------------------------------------------

# The livecyte dataset (livecyte_data.tsv) contains per-frame
# measurements for each automatically tracked cell, with columns
# including clone, replicate, tracking.id, frame, position.x,
# position.y, track.length, displacement, instantaneous.velocity,
# and morphological features (volume, radius, dry.mass, etc.).
#
# The manual dataset (manual_data.tsv) contains one row per manually
# tracked cell, with columns: LID, TID, track.duration, track.length,
# meandering.index, euclidean.distance, mean.speed, and cell.line
# (A1 or B2).


# Analysis overview ----------------------------------------------

# This script collapses the per-frame livecyte data to one row per
# tracking.id, then matches each manual TID to the closest livecyte
# tracking.id using three approximately-equal metrics:
#    track.length       ~= total_path_length
#    euclidean.distance  ~= final_displacement
#    mean.speed (um/min) ~= mean.speed * 60  (livecyte is um/sec)
#
# The distance function is a weighted sum of absolute ratio errors
# with weights 2.0 for path length, 2.0 for displacement, and 0.1
# for speed. Speed gets lower weight because it diverges more between
# the two tracking methods. The output is a CSV file that maps each
# TID to its corresponding tracking.id.

# ================================================================

# Packages required ----------------------------------------------

# for data manipulation and plotting
library(tidyverse)


# Data import ----------------------------------------------------

# livecyte per-frame dataset
livecyte <- read_tsv("project/data/movement_morphology/livecyte_data.tsv",
                     show_col_types = FALSE)

# manual tracking dataset
manual <- read_tsv("project/data/movement_morphology/manual_data.tsv",
                   show_col_types = FALSE)


# Collapse livecyte to one row per tracking.id -------------------

# Each tracking.id has many per-frame rows. We collapse to summary
# metrics so that each tracking.id has a single total_path_length,
# final_displacement, and mean.speed to compare against the manual
# data
livecyte_collapsed_unfiltered <- livecyte |>
  arrange(clone, replicate, tracking.id, frame) |>
  group_by(clone, replicate, tracking.id) |>
  summarise(
    n_frames = n(),
    start_frame = first(frame),
    total_path_length = last(track.length),
    final_displacement = last(displacement),
    volume = mean(volume),
    radius = mean(radius),
    mean.thickness = mean(mean.thickness),
    sphericity = mean(sphericity),
    length.to.width.ratio = mean(length.to.width),
    dry.mass = mean(dry.mass),
    mean.speed = mean(instantaneous.velocity, na.rm = TRUE) * 60,
    .groups = "drop"
  )
# mean.speed is converted from um/sec to um/min by multiplying by 60

# Remove NaN mean.speed entries. These are likely debris, as they
# have no speed (e.g. tracked for only one frame)
livecyte_collapsed_unfiltered <- livecyte_collapsed_unfiltered |>
  filter(!is.na(mean.speed))

# Save this collapsed dataset for use in later scripts
write_tsv(livecyte_collapsed_unfiltered,
          "project/data/movement_morphology/livecyte_collapsed_unfiltered.tsv")


# Filter to replicates present in manual data --------------------

# The manual dataset only contains cell.line A1 (clone A replicate 1)
# and B2 (clone B replicate 2). We remove the other four
# clone-replicate combinations
livecyte_A1_B2 <- livecyte_collapsed_unfiltered |>
  filter(
    (clone == "A" & replicate == 1) |
      (clone == "B" & replicate == 2)
  ) |>
  mutate(
    cell.line = paste0(clone, replicate),
    mean.speed.converted = mean.speed * 60
  )
# mean.speed.converted is now in um/hr for the matching step;
# manual mean.speed is in um/min


# Distance function ----------------------------------------------

# We use the absolute ratio error: |a - b| / max(a, b). This is
# bounded between 0 and 1, and handles cases where one value is much
# larger than the other without being dominated by scale differences
ratio_err <- function(a, b) {
  abs(a - b) / pmax(a, b, 1e-6)
}


# Matching function ----------------------------------------------

# For each manual track, we score every candidate livecyte row and
# pick the one with the lowest weighted error. Path length and
# displacement get weight 2.0 each as they are the most reliable
# metrics; speed gets 0.1 as it is noisier between the two methods
w_length <- 2.0
w_disp   <- 2.0
w_speed  <- 0.1

match_one <- function(m_row, candidates) {
  m_len  <- m_row$track.length
  m_disp <- m_row$euclidean.distance
  m_spd  <- m_row$mean.speed
  
  scored <- candidates |>
    mutate(
      err_length = ratio_err(total_path_length, m_len),
      err_disp   = ratio_err(final_displacement, m_disp),
      err_speed  = ratio_err(mean.speed.converted, m_spd),
      total_err  = w_length * err_length +
        w_disp * err_disp +
        w_speed * err_speed
    )
  
  best <- scored |>
    slice_min(total_err, n = 1, with_ties = FALSE)
  
  tibble(
    tracking.id = best$tracking.id,
    err_length  = best$err_length,
    err_disp    = best$err_disp,
    err_speed   = best$err_speed,
    total_err   = best$total_err
  )
}

# Run matching across every manual track
results <- manual |>
  rowwise() |>
  mutate(
    match = list(
      match_one(
        cur_data_all(),
        livecyte_A1_B2 |> filter(cell.line == cur_data_all()$cell.line)
      )
    )
  ) |>
  unnest(match) |>
  ungroup()


# Save result ----------------------------------------------------

# Build final output with LID, TID, clone, replicate, tracking.id,
# and cell.line
mapping <- results |>
  mutate(
    clone     = str_sub(cell.line, 1, 1),
    replicate = as.integer(str_sub(cell.line, 2, 2))
  ) |>
  select(LID, TID, clone, replicate, tracking.id, cell.line)

write_csv(mapping,
          "project/data/movement_morphology/TID_to_trackingid.csv")


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