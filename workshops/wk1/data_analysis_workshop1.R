# Data analysis workshop 1

# Packages -------------------------------
library(tidyverse)    # For data manipulation

# Load data -------------------------------
data <- read_tsv("workshops/wk1/practice_dataset.tsv")

# Data manipulation -------------------------------
cloneA <- data |>                 # Create a new data frame for cloneA
  filter(clone == "cloneA") |> 
  select(clone,                   # Select only the columns we need
         replicate,
         length.to.width)

cloneB <- data |>                 # Create a new data frame for cloneB 
  filter(clone == "cloneB") |> 
  select(clone,                   # Select only the columns we need
         replicate,
         length.to.width)

cloneA |> 
  group_by(replicate) |>           # Group the data by replicate
  summarise(mean_length_to_width = mean(length.to.width)) # Calculate the mean length to width ratio for each replicate

cloneB |>
  group_by(replicate) |>           # Group the data by replicate
  summarise(mean_length_to_width = mean(length.to.width)) # Calculate the mean length to width ratio for each replicate

# Statistics -------------------------------
# Test whether ltw is normally distributed by using a ks test
ks.test(cloneA$length.to.width, "pnorm", mean = mean(cloneA$length.to.width), sd = sd(cloneA$length.to.width))
ks.test(cloneB$length.to.width, "pnorm", mean = mean(cloneB$length.to.width), sd = sd(cloneB$length.to.width))

# p-values for both are less than 0.05, so we can reject the null hypothesis that the data is normally distributed
# Therefore we will use a non-parametric test to compare the two clones
wilcox.test(cloneA$length.to.width, cloneB$length.to.width)
# The p-value is less than 0.05, so we can reject the null hypothesis that the two clones have the same length to width ratio
# Therefore there is evidence to suggest that the two clones have different length to width ratios

# Visualisation -------------------------------
ggplot(data, aes(x = clone, y = length.to.width)) +
  geom_boxplot() +
  geom_violin() +
  labs(title = "Length to Width Ratio of Clones A and B",
       x = "Clone",
       y = "Length to Width Ratio") +
  theme_minimal()
