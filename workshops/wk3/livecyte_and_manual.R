# Workshop 3: Automated and manual livecyte data analysis

# Packages
library(tidyverse)
library(ggpubr)
library(corrr)

# Load data
livecyte <-read_tsv(url("https://djeffares.github.io/BIO66I/data/all-cell-data-FFT.filtered.2024-02-22.tsv"),
                 col_types = cols(
                   clone = col_factor(),
                   replicate = col_factor(),
                   tracking.id=col_factor(),
                   lineage.id=col_factor()
                 )
)

# Select only the columns we need
cell.movement.data <- select(livecyte,
                             clone,
                             replicate,
                             displacement, 
                             track.length, 
                             instantaneous.velocity
)


# Now we will look at summarising the data

# Plot the data
ggplot(cell.movement.data,aes(x = clone,y = log10(displacement),colour = replicate))+
  geom_violin(alpha=0.5)+
  stat_compare_means()

# Split data by clone
cloneA.data <- cell.movement.data |> filter(clone == "cloneA")
cloneB.data <- cell.movement.data |> filter(clone == "cloneB")

numeric.columns <- cell.movement.data |> 
  select(where(is.numeric)) |> 
  names()

# Create an empty dataframe to deposit our summarised results into
clone.comparisons <- data.frame(
  variable = character(),
  cloneA.median = numeric(),
  cloneB.median = numeric(),
  median.ratio = numeric(),
  p.value = numeric()
)

# Loop through the variables we want to compare
for(column.name in numeric.columns) {
  
  #calculate median for clone A
  cloneA.median <- median(cloneA.data[[column.name]], na.rm = TRUE)
  
  #calculate median for clone B  
  cloneB.median <- median(cloneB.data[[column.name]], na.rm = TRUE)
  
  #calculate the median ratio (cloneA.median / cloneB.median)
  ratio <- cloneA.median / cloneB.median
  
  #run wilcox test comparing the two clones for this variable
  test.result <- wilcox.test(livecyte[[column.name]] ~ livecyte$clone)
  
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

