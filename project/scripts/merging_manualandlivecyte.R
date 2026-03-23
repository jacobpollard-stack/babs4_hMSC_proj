#!/usr/bin/env Rscript

# Load required libraries
library(tidyverse)
library(data.table)

# ============================================================================
# 1. READ AND PREPARE DATA
# ============================================================================

cat("Loading datasets...\n")

# Read livecyte data (already has clone and replicate columns)
livecyte <- read.csv('/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 2/BABS/babs4_hMSC_proj/project/data/movement_morphology/livecyte_pretty2.csv', row.names = 1)

# Read manual data (has cell.line column that needs conversion)
manual <- read.delim('/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 2/BABS/babs4_hMSC_proj/project/data/movement_morphology/manual_data.tsv')

cat("✓ Livecyte data loaded:", nrow(livecyte), "rows x", ncol(livecyte), "cols\n")
cat("✓ Manual data loaded:", nrow(manual), "rows x", ncol(manual), "cols\n")

# ============================================================================
# 2. CONVERT CELL LINE TO CLONE/REPLICATE IN MANUAL DATA
# ============================================================================

# Parse cell.line (A1, B2) into clone and replicate
manual <- manual %>%
  mutate(
    clone = case_when(
      substr(cell.line, 1, 1) == "A" ~ "cloneA",
      substr(cell.line, 1, 1) == "B" ~ "cloneB",
      TRUE ~ NA_character_
    ),
    replicate = as.numeric(substr(cell.line, 2, 3))
  )

cat("\nManual data contains:\n")
manual_groups <- manual %>% distinct(cell.line, clone, replicate)
print(as.data.frame(manual_groups))

cat("\nLivecyte data contains:\n")
livecyte_groups <- livecyte %>% distinct(clone, replicate) %>% arrange(clone, replicate)
print(as.data.frame(livecyte_groups))

# ============================================================================
# 3. NORMALIZE PARAMETERS FOR MATCHING
# ============================================================================

# The three parameters to match:
# livecyte: total_path_length, final_displacement, mean.speed
# manual: track.length, euclidean.distance, mean.speed

# Function to normalize values to 0-1 scale for fair distance calculation
normalize <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  if (min_val == max_val) {
    return(rep(0.5, length(x)))  # Handle case where all values are identical
  }
  (x - min_val) / (max_val - min_val)
}

cat("\nNormalizing parameters within each clone/replicate group...\n")

livecyte <- livecyte %>%
  group_by(clone, replicate) %>%
  mutate(
    total_path_length_norm = normalize(total_path_length),
    final_displacement_norm = normalize(final_displacement),
    mean.speed_norm = normalize(mean.speed)
  ) %>%
  ungroup()

manual <- manual %>%
  group_by(clone, replicate) %>%
  mutate(
    track.length_norm = normalize(track.length),
    euclidean.distance_norm = normalize(euclidean.distance),
    mean.speed_norm = normalize(mean.speed)
  ) %>%
  ungroup()

# ============================================================================
# 4. MATCHING FUNCTION
# ============================================================================

# For each cell line group, find best matches using 3D Euclidean distance
match_objects <- function(livecyte_subset, manual_subset) {
  
  matches_list <- list()
  
  for (i in seq_len(nrow(livecyte_subset))) {
    
    live_row <- livecyte_subset[i, ]
    
    # Calculate 3D Euclidean distance to all manual rows in this group
    distances <- sqrt(
      (manual_subset$track.length_norm - live_row$total_path_length_norm)^2 +
        (manual_subset$euclidean.distance_norm - live_row$final_displacement_norm)^2 +
        (manual_subset$mean.speed_norm - live_row$mean.speed_norm)^2
    )
    
    # Create results dataframe for this livecyte object
    match_results <- data.frame(
      livecyte_tracking_id = live_row$tracking.id,
      livecyte_lineage_id = live_row$lineage.id,
      manual_tracking_id = manual_subset$TID,
      manual_lineage_id = manual_subset$LID,
      distance = distances,
      rank = rank(distances),
      stringsAsFactors = FALSE
    )
    
    matches_list[[i]] <- match_results
  }
  
  return(bind_rows(matches_list))
}

# ============================================================================
# 5. PERFORM MATCHING FOR EACH CELL LINE GROUP
# ============================================================================

cat("\nPerforming matching for A1 and B2 only...\n")

all_matches <- data.frame()

# Only match the cell lines present in manual data (A1 and B2)
target_groups <- manual %>% 
  distinct(clone, replicate) %>% 
  arrange(clone, replicate)

for (i in seq_len(nrow(target_groups))) {
  
  clone_val <- target_groups$clone[i]
  rep_val <- target_groups$replicate[i]
  
  live_subset <- livecyte %>% 
    filter(clone == clone_val & replicate == rep_val)
  
  man_subset <- manual %>% 
    filter(clone == clone_val & replicate == rep_val)
  
  if (nrow(live_subset) > 0 & nrow(man_subset) > 0) {
    
    cat("  Matching", clone_val, "replicate", rep_val, ":", 
        nrow(live_subset), "livecyte objects to", 
        nrow(man_subset), "manual objects\n")
    
    matches <- match_objects(live_subset, man_subset)
    matches$clone <- clone_val
    matches$replicate <- rep_val
    
    all_matches <- bind_rows(all_matches, matches)
  } else {
    cat("  ⚠ Skipping", clone_val, "replicate", rep_val, 
        "(missing from one or both datasets)\n")
  }
}

# ============================================================================
# 6. IDENTIFY BEST MATCHES AND FLAGS
# ============================================================================

cat("\nEvaluating match quality...\n")

best_matches <- all_matches %>%
  group_by(livecyte_tracking_id, livecyte_lineage_id) %>%
  arrange(distance) %>%
  mutate(
    best_match = (rank == 1),
    match_quality = case_when(
      distance < 0.3 ~ "Excellent",
      distance < 0.5 ~ "Good",
      distance < 0.7 ~ "Acceptable",
      TRUE ~ "Poor"
    ),
    num_close_matches = sum(distance < 0.5),
    num_very_close = sum(distance < 0.3)
  ) %>%
  ungroup()

# ============================================================================
# 7. CREATE MERGED DATASET
# ============================================================================

cat("\nCreating merged dataset...\n")

# Get best matches only
best_only <- best_matches %>%
  filter(best_match == TRUE) %>%
  select(-best_match, -rank)

# Merge with original livecyte data
merged_dataset <- best_only %>%
  left_join(
    livecyte %>% select(clone, replicate, tracking.id, lineage.id, everything()),
    by = c("clone" = "clone", 
           "replicate" = "replicate", 
           "livecyte_tracking_id" = "tracking.id",
           "livecyte_lineage_id" = "lineage.id")
  ) %>%
  left_join(
    manual %>% select(clone, replicate, TID, LID, everything()),
    by = c("clone" = "clone", 
           "replicate" = "replicate", 
           "manual_tracking_id" = "TID",
           "manual_lineage_id" = "LID"),
    suffix = c("_livecyte", "_manual")
  )

# ============================================================================
# 8. SUMMARY STATISTICS AND OUTPUT
# ============================================================================

cat("\n")
cat("=" %s>% rep("=", 70) %s>% paste(collapse = ""), "\n")
cat("MATCHING SUMMARY (cloneA/replicate 1 and cloneB/replicate 2)\n")
cat("=" %s>% rep("=", 70) %s>% paste(collapse = ""), "\n\n")

cat("Livecyte dataset summary:\n")
cat("  - cloneA: replicates 1, 2, 3 (", 
    nrow(livecyte %>% filter(clone == "cloneA")), " total objects)\n")
cat("  - cloneB: replicates 1, 2, 3 (", 
    nrow(livecyte %>% filter(clone == "cloneB")), " total objects)\n\n")

cat("Manual dataset summary (matching targets):\n")
cat("  - cloneA/replicate 1: ", nrow(manual %>% filter(clone == "cloneA" & replicate == 1)), " objects\n")
cat("  - cloneB/replicate 2: ", nrow(manual %>% filter(clone == "cloneB" & replicate == 2)), " objects\n")
cat("  - Total manual objects: ", nrow(manual), "\n\n")

cat("Matching results:\n")
cat("  - cloneA/replicate 1 livecyte objects: ", 
    nrow(livecyte %>% filter(clone == "cloneA" & replicate == 1)), "\n")
cat("  - cloneB/replicate 2 livecyte objects: ", 
    nrow(livecyte %>% filter(clone == "cloneB" & replicate == 2)), "\n")
cat("  - Successful matches: ", nrow(best_matches %>% filter(best_match == TRUE)), "\n")
cat("  - Match rate: ", 
    round(nrow(best_matches %>% filter(best_match == TRUE)) / nrow(livecyte %>% filter((clone == "cloneA" & replicate == 1) | (clone == "cloneB" & replicate == 2))) * 100, 1), 
    "%\n\n")

# Quality breakdown
cat("Match Quality Distribution:\n")
quality_dist <- best_matches %>% 
  filter(best_match == TRUE) %>% 
  count(match_quality) %>%
  arrange(factor(match_quality, levels = c("Excellent", "Good", "Acceptable", "Poor")))
print(as.data.frame(quality_dist))

# Objects with multiple close candidates
multiples <- best_matches %>% 
  filter(best_match == TRUE & num_close_matches > 1)

cat("\nObjects with multiple close candidates (<0.5 distance):", nrow(multiples), "\n")
if (nrow(multiples) > 0) {
  cat("\nTop 20 objects with multiple candidates:\n")
  print(multiples %>% 
          select(clone, replicate, livecyte_tracking_id, manual_tracking_id, 
                 num_close_matches, num_very_close, distance) %>%
          arrange(desc(num_close_matches), distance) %>%
          head(20) %>%
          as.data.frame())
}

# ============================================================================
# 9. OUTPUT FILES
# ============================================================================

cat("\n")
cat("=" %s>% rep("=", 70) %s>% paste(collapse = ""), "\n")
cat("WRITING OUTPUT FILES\n")
cat("=" %s>% rep("=", 70) %s>% paste(collapse = ""), "\n\n")

output_dir <- '/Users/jacobpollard/Documents/Uni/Biology/Second year/Sem 2/BABS/babs4_hMSC_proj/project'
dir.create(output_dir, showWarnings = FALSE)

# Write merged dataset
output_path_merged <- file.path(output_dir, "merged_dataset.csv")
write.csv(merged_dataset, output_path_merged, row.names = FALSE)
cat("✓ Merged dataset saved to: merged_dataset.csv\n")
cat("  Dimensions:", nrow(merged_dataset), "rows x", ncol(merged_dataset), "cols\n")

# Write detailed match report (all candidates)
output_path_all <- file.path(output_dir, "match_report_all_candidates.csv")
write.csv(best_matches, output_path_all, row.names = FALSE)
cat("✓ Detailed match report saved to: match_report_all_candidates.csv\n")
cat("  Dimensions:", nrow(best_matches), "rows x", ncol(best_matches), "cols\n")

# Write best matches summary
output_path_best <- file.path(output_dir, "best_matches_only.csv")
write.csv(best_only, output_path_best, row.names = FALSE)
cat("✓ Best matches summary saved to: best_matches_only.csv\n")
cat("  Dimensions:", nrow(best_only), "rows x", ncol(best_only), "cols\n")

# ============================================================================
# 10. VISUAL SUMMARY
# ============================================================================

cat("\n")
cat("=" %s>% rep("=", 70) %s>% paste(collapse = ""), "\n")
cat("TOP 15 BEST MATCHES (LOWEST DISTANCE)\n")
cat("=" %s>% rep("=", 70) %s>% paste(collapse = ""), "\n\n")

top_matches <- all_matches %>%
  arrange(distance) %>%
  head(15) %>%
  left_join(best_matches %>% filter(best_match == TRUE) %>% select(livecyte_tracking_id, match_quality),
            by = "livecyte_tracking_id")

print(
  top_matches %>%
    select(clone, replicate, livecyte_tracking_id, manual_tracking_id, 
           distance, rank, match_quality) %>%
    as.data.frame()
)

cat("\n")
cat("=" %s>% rep("=", 70) %s>% paste(collapse = ""), "\n")
cat("WORST 15 BEST MATCHES (HIGHEST DISTANCE)\n")
cat("=" %s>% rep("=", 70) %s>% paste(collapse = ""), "\n\n")

worst_matches <- best_matches %>%
  filter(best_match == TRUE) %>%
  arrange(desc(distance)) %>%
  head(15)

print(
  worst_matches %>%
    select(clone, replicate, livecyte_tracking_id, manual_tracking_id, 
           distance, match_quality, num_close_matches) %>%
    as.data.frame()
)

cat("\n")
cat("=" %s>% rep("=", 70) %s>% paste(collapse = ""), "\n")
cat("SCRIPT COMPLETE\n")
cat("=" %s>% rep("=", 70) %s>% paste(collapse = ""), "\n")

