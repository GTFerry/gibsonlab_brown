library(dplyr)
library(tibble)
library(tidyr)

# Load data
data <- read.csv("results/average_data.csv", row.names = 1)

# First convert rownames to a column named cluster, then extract the values of each param and create a column for each.
data <- data %>%
  rownames_to_column(var = "cluster") %>%
  mutate(
    ResRNA = gsub(".*ResRNA_([0-9.]+).*", "\\1", cluster),
    ResATAC = gsub(".*ResATAC_([0-9.]+).*", "\\1", cluster),
    Iters = gsub(".*Iters_([0-9]+).*", "\\1", cluster),
    Dims = gsub(".*Dims_(1to[0-9]+).*", "\\1", cluster)
  )

# reshape data from wide format to long format
mean_scores <- data %>%
  pivot_longer(cols = c(ResRNA, ResATAC, Iters, Dims), names_to = "parameter", values_to = "value") %>%
  group_by(parameter, value) %>%
  summarize(mean_score = mean(x), .groups = 'drop')

best_parameters <- mean_scores %>%
  group_by(parameter) %>%
  slice(which.max(mean_score))

print(best_parameters)
