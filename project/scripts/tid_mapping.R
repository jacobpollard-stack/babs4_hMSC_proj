# ====================================================================
# Livecyte Data Analysis
# - Mapping manual dataset to the livecyte dataset
# Please turn on Soft Wrap in your code editor to view the full code without horizontal scrolling.
#
# In the experiment, Livecyte used its automatic cell tracking software to assign a variety of metrics to each cell track. This tracking is very noisy, and the tracking IDs (tracking.id) do not directly correspond to the manual tracking IDs (TID) in the manual ImageJ dataset. In order to perform a direct comparison between the two datasets, we need to match each TID to the corresponding tracking.id. This code performs this by comparing the track length and final displacement between the two datasets, and finding the closest match for each TID. The output is a CSV file that maps each TID to its corresponding tracking.id.
#
# ====================================================================
#
# 1. Load libraries --------------------------------------------------
#
library(tidyverse)
#
# 2. Read data -------------------------------------------------------
#
# 2a. Load the livecyte dataset
#
livecyte <- read_tsv("project/data/movement_morphology/livecyte_collapsed.tsv", show_col_types = FALSE)
#
# 2b. Load the manual dataset
#
manual   <- read_tsv("project/data/movement_morphology/manual_data.tsv", show_col_types = FALSE)

# 2. Filter livecyte to keep only the replicates in the manually tracked dataset ------------------------------------------------------
#
# Remove: clone A replicate 2, clone A replicate 3,
#         clone B replicate 1, clone B replicate 3
#
livecyte_A1_B2 <- livecyte |> 
  filter(
    (clone == "A" & replicate == 1) |
      (clone == "B" & replicate == 2)
  ) |> 
  mutate(cell.line = paste0(clone, replicate))
#
# 3. Define distance function ------------------------
#
# We will use the squared relative error on track.length and euclidean.distance
# mean.speed is too noisy between datasets to use
#
# This function computes the squared relative error, with a small epsilon to avoid dividing by zero. We square the error to penalise large errors more heavily
#
compute_rel_err <- function(manual_val, livecyte_val) {
  ((livecyte_val - manual_val) / pmax(abs(manual_val), 1e-6))^2
}
#
# 4. Building the matching function ------------------
#
match_one <- function(m_row, candidates) {
  m_len  <- m_row$track.length
  m_disp <- m_row$euclidean.distance
  m_spd  <- m_row$mean.speed
  
  scored <- candidates |> 
    mutate(
      err_length = compute_rel_err(m_len,  total_path_length),
      err_disp   = compute_rel_err(m_disp, final_displacement),
      err_speed  = compute_rel_err(m_spd,  mean.speed),
      total_err  = err_length + err_disp # length + displacement only
    )

  best <- scored |> 
    slice_min(total_err, n = 1, with_ties = FALSE)
  
  tibble(
    tracking.id = best$tracking.id,
    err_length  = sqrt(best$err_length),
    err_disp    = sqrt(best$err_disp),
    err_speed   = sqrt(best$err_speed),
    total_err   = best$total_err
  )
}
#
# Match for each manual track
#
results <- manual |> 
  rowwise() |> 
  mutate(
    match = list(
      match_one(
        cur_data_all(),
        livecyte_A1_B2 |>  filter(cell.line == cur_data_all()$cell.line)
      )
    )
  ) |> 
  unnest(match) |> 
  ungroup()
#
# 5. Save result -------------------------------------
#
# 5a. Create a final output with LID, TID, clone, replicate, tracking.id, and cell.line
#
final_output <- results |> 
  mutate(
    clone     = str_sub(cell.line, 1, 1),
    replicate = as.integer(str_sub(cell.line, 2, 2))
  ) |> 
  select(LID, TID, clone, replicate, tracking.id, cell.line)
#
# 5b. Save the output as a CSV file
#
write_csv(final_output, "project/data/movement_morphology/TID_to_trackingid.csv")
