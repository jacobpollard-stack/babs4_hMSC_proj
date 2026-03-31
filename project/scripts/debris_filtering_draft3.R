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
library(kSamples)
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
# Clone A appears more likely for the automated Livecyte tracking to lose cells, leading to a different distributions between manual and livecyte for A1 than B2.
#
# 3c. Perform Wilcoxon rank sum and Kolmogorov-Smirnov tests to determine if the distributions of the manual and livecyte_filtered datasets are significantly different for each parameter. This will allow us to quantify the accuracy of the livecyte tracking data.
# A1 euclidean distance vs final displacement
w_A1_euc <- wilcox.test(manual |> filter(cell.line == "A1") |> pull(euclidean.distance),
                        livecyte_filtered |> filter(clone == "cloneA", replicate == "1") |> pull(final_displacement))
k_A1_euc <- ks.test(manual |> filter(cell.line == "A1") |> pull(euclidean.distance),
                    livecyte_filtered |> filter(clone == "cloneA", replicate == "1") |> pull(final_displacement))
#
# A1 track length vs total path length
w_A1_trk <- wilcox.test(manual |> filter(cell.line == "A1") |> pull(track.length),
                        livecyte_filtered |> filter(clone == "cloneA", replicate == "1") |> pull(total_path_length))
k_A1_trk <- ks.test(manual |> filter(cell.line == "A1") |> pull(track.length),
                    livecyte_filtered |> filter(clone == "cloneA", replicate == "1") |> pull(total_path_length))
#
# A1 mean speed vs mean speed
w_A1_spd <- wilcox.test(manual |> filter(cell.line == "A1") |> pull(mean.speed),
                        livecyte_filtered |> filter(clone == "cloneA", replicate == "1") |> pull(mean.speed))
k_A1_spd <- ks.test(manual |> filter(cell.line == "A1") |> pull(mean.speed),
                    livecyte_filtered |> filter(clone == "cloneA", replicate == "1") |> pull(mean.speed))
#
# B2 euclidean distance vs final displacement
w_B2_euc <- wilcox.test(manual |> filter(cell.line == "B2") |> pull(euclidean.distance),
                        livecyte_filtered |> filter(clone == "cloneB", replicate == "2") |> pull(final_displacement))
k_B2_euc <- ks.test(manual |> filter(cell.line == "B2") |> pull(euclidean.distance),
                    livecyte_filtered |> filter(clone == "cloneB", replicate == "2") |> pull(final_displacement))
#
# B2 track length vs total path length
w_B2_trk <- wilcox.test(manual |> filter(cell.line == "B2") |> pull(track.length),
                        livecyte_filtered |> filter(clone == "cloneB", replicate == "2") |> pull(total_path_length))
k_B2_trk <- ks.test(manual |> filter(cell.line == "B2") |> pull(track.length),
                    livecyte_filtered |> filter(clone == "cloneB", replicate == "2") |> pull(total_path_length))
#
# B2 mean speed vs mean speed
w_B2_spd <- wilcox.test(manual |> filter(cell.line == "B2") |> pull(mean.speed),
                        livecyte_filtered |> filter(clone == "cloneB", replicate == "2") |> pull(mean.speed))
k_B2_spd <- ks.test(manual |> filter(cell.line == "B2") |> pull(mean.speed),
                    livecyte_filtered |> filter(clone == "cloneB", replicate == "2") |> pull(mean.speed))
#
# Collect into tibble
test_results <- tibble(
  cell_line  = c("A1",      "A1",       "A1",      "B2",      "B2",       "B2"),
  parameter  = c("displacement", "path_length", "mean_speed",
                 "displacement", "path_length", "mean_speed"),
  wilcox_W   = c(w_A1_euc$statistic, w_A1_trk$statistic, w_A1_spd$statistic,
                 w_B2_euc$statistic, w_B2_trk$statistic, w_B2_spd$statistic),
  wilcox_p   = c(w_A1_euc$p.value,   w_A1_trk$p.value,   w_A1_spd$p.value,
                 w_B2_euc$p.value,   w_B2_trk$p.value,   w_B2_spd$p.value),
  ks_D       = c(k_A1_euc$statistic, k_A1_trk$statistic, k_A1_spd$statistic,
                 k_B2_euc$statistic, k_B2_trk$statistic, k_B2_spd$statistic),
  ks_p       = c(k_A1_euc$p.value,   k_A1_trk$p.value,   k_A1_spd$p.value,
                 k_B2_euc$p.value,   k_B2_trk$p.value,   k_B2_spd$p.value)
)
#
# Displacement is the only reliable debris-filtering parameter. The distributions for the other two parameters vary too much between the manual and livecyte datasets, likely due to the livecyte tracking losing cells, especially in A1.
#
# 4. Check whether the debris is more likely to be tracked in certain clones or replicates by plotting the distribution of final_displacement for each replicate in the livecyte_pretty dataset.
#
# 4a. Split datasets into A1 and B2
pretty_A <- livecyte_pretty |> filter(clone == "cloneA" & replicate == "1")
pretty_B <- livecyte_pretty |> filter(clone == "cloneB" & replicate == "2")
#
# 4b. Create a 'debris' dataset from livecyte_pretty using livecyte_filtered
debris_A <- pretty_A |> 
anti_join(tid_mapping, by = c("clone", "replicate", "tracking.id"))
debris_B <- pretty_B |> 
  anti_join(tid_mapping, by = c("clone", "replicate", "tracking.id"))
#
# 4c. Plot the distribution of the debris in A1 and B2
# A
ggplot(debris_A, aes(x = final_displacement)) +
  geom_histogram(binwidth = 10)
# B
ggplot(debris_B, aes(x = final_displacement)) +
  geom_histogram(binwidth = 10)
#
# 4d. Check homogeneity of debris across replicates using statistical tests.
# 4di. Create full debris dataset across all replicates
livecyte_debris <- livecyte_pretty |>
  anti_join(tid_mapping, by = c("clone", "replicate", "tracking.id"))
#
# 4dii. Split by clone and drop unused factor levels
debris_A_all <- livecyte_debris |> filter(clone == "cloneA") |> droplevels()
debris_B_all <- livecyte_debris |> filter(clone == "cloneB") |> droplevels()
#
# 4diii. Perform AD test to check for homogeneity of debris distribution across replicates within each clone
#
# A
ad_debris_A <- do.call(ad.test, split(debris_A_all$final_displacement, debris_A_all$replicate))
#
# B
ad_debris_B <- do.call(ad.test, split(debris_B_all$final_displacement, debris_B_all$replicate))
#
# 4div. Perform pairwise Wilcoxon tests for pairwise differences, with BH correction
pairwise.wilcox.test(debris_A_all$final_displacement, debris_A_all$replicate, p.adjust.method = "BH")
pairwise.wilcox.test(debris_B_all$final_displacement, debris_B_all$replicate, p.adjust.method = "BH")
#
# Clone B replicate 3 seems to be highly significantly different from the other replicates, which could be due to a biological difference, or a difference in debris distribution.
#
#  4e. Plot the distribution of the debris in each replicate to visually check for differences.
# A
ggplot(debris_A_all, aes(x = final_displacement)) +
  geom_histogram(binwidth = 10) +
  facet_wrap(~ replicate)
# B
ggplot(debris_B_all, aes(x = final_displacement)) +
  geom_histogram(binwidth = 10) +
  facet_wrap(~ replicate)
#
# This is too much work for too little gain. The debris distribution is not homogenous across replicates, but it's difficult to determine how this overlaps with cells.