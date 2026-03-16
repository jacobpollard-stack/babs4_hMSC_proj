# Workshop 2: working with the real dataset

# Packages -------------------------------------
library(tidyverse)
library(ggpubr)

# Load the data -------------------------------------
cells <-read_tsv(url("https://djeffares.github.io/BIO66I/data/all-cell-data-FFT.filtered.2024-02-22.tsv"),
                 col_types = cols(
                   clone = col_factor(),
                   replicate = col_factor(),
                   tracking.id = col_factor(),
                   lineage.id = col_factor()
                 )
)

# Remove all the movement metrics (which are not reliable)
cells <- select(cells, 
                -position.x,
                -position.y, 
                -pixel.position.x, 
                -pixel.position.y,
                -instantaneous.velocity,
                -instantaneous.velocity.x,
                -instantaneous.velocity.y,
                -track.length
)

# Check that our data is simpler
names(cells)

str(cells)

# Plot the data -------------------------------------
ggplot(cells, aes(x = clone, y = displacement, fill = replicate)) +
  geom_boxplot() +
  theme_bw() +
  geom_violin()


# Make cloneA and cloneB data frames
cloneA.data <- cells |> filter(clone == "cloneA")
cloneB.data <- cells |> filter(clone == "cloneB")

#get the names of the numeric columns
numeric.columns <- cells |> 
  select(where(is.numeric)) |> 
  names()

#see what we have
numeric.columns

#create empty data frame
clone.comparisons <- data.frame(
  variable = character(),
  cloneA.median = numeric(),
  cloneB.median = numeric(),
  median.ratio = numeric(),
  p.value = numeric()
)

#loop through each numeric column
for(column.name in numeric.columns) {
  
  #calculate median for clone A
  cloneA.median <- median(cloneA.data[[column.name]], na.rm = TRUE)
  
  #calculate median for clone B  
  cloneB.median <- median(cloneB.data[[column.name]], na.rm = TRUE)
  
  #calculate the median ratio (cloneA.median / cloneB.median)
  ratio <- cloneA.median / cloneB.median
  
  #run wilcox test comparing the two clones for this variable
  test.result <- wilcox.test(cells[[column.name]] ~ cells$clone)
  
  #extract the p-value
  p.val <- test.result$p.value
  
  #add this row to our results table
  new.row <- data.frame(
    variable = column.name,
    cloneA.median = signif(cloneA.median,2),
    cloneB.median = signif(cloneB.median,2),
    median.ratio = signif(ratio,2),
    p.value = p.val
  )
  
  #add the new.row to the clone.comparisons data frame
  clone.comparisons <- rbind(clone.comparisons, new.row)
}


# Principal component analysis -------------------------------------
# Select only the numeric columns for PCA
cells_for_pca <- cells |>
  #select the numeric columns
  select(clone, replicate, volume, mean.thickness, radius, area, 
         sphericity, length, width, orientation, dry.mass, length.to.width) |>
  # remove rows with NA values
  drop_na()

pca_result <- cells_for_pca |>
  select(-clone, -replicate) |>  # remove the group categories
  prcomp()                        # performs the PCA

# View summary of PCA
summary(pca_result)

pca_scores <- as.data.frame(pca_result$x) |>
  bind_cols(cells_for_pca |> select(clone, replicate))

pca.plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = clone, fill = clone)) +
  geom_point(alpha = 0.1, size = 4) 

# View the plot
pca.plot

pca.plot +
  stat_ellipse(geom = "polygon", alpha = 0.1, linewidth = 1.2) +
  scale_color_manual(values = c("cloneA" = "red", "cloneB" = "blue"))+
  facet_wrap(~replicate)+
  theme_classic2()


library(readr)

url <- "https://djeffares.github.io/BIO66I/data/A1-and-B2-tracking.data.2025-02-27.tsv"
data <- read_tsv(url)

write_tsv(data, "all-cell-data-FFT.filtered.2024-02-22.tsv")
