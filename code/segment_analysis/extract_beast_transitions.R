#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  # Default paths if not provided
  tree_file <- "results/beast/final_run/H5N1_HA_panel_postQC.mcc.tree"
  out_file <- "results/beast/final_run/joint_transitions_summary.csv"
} else {
  tree_file <- args[1]
  out_file <- args[2]
}

log_file <- "results/beast/final_run/H5N1_HA_panel_postQC.log"

suppressPackageStartupMessages({
  library(treeio)
  library(dplyr)
  library(tidyr)
})

cat("Reading BEAST log to calculate Bayes Factors...\n")
log_data <- read.table(log_file, header=TRUE, comment.char="#", sep="\t", check.names=FALSE)
burnin <- floor(0.1 * nrow(log_data))
log_data <- log_data[-(1:burnin), ]

# Helper to find BSSVS indicator column safely (since they are ordered alphabetically by BEAST)
get_indicator_mean <- function(trait, state1, state2) {
  col1 <- paste0(trait, ".indicators.", state1, ".", state2)
  col2 <- paste0(trait, ".indicators.", state2, ".", state1)
  
  if (col1 %in% colnames(log_data)) return(mean(log_data[[col1]]))
  if (col2 %in% colnames(log_data)) return(mean(log_data[[col2]]))
  return(NA)
}

# Calculate BFs
# Location: K=13, E=12, total=78. q = 12/78
q_loc <- 12 / 78
prior_odds_loc <- q_loc / (1 - q_loc)

# Host: K=3, E=2, total=3. q = 2/3
q_host <- 2 / 3
prior_odds_host <- q_host / (1 - q_host)

cat("Reading tree...\n")
tree <- read.beast(tree_file)
df <- as_tibble(tree)

lookup <- df %>% select(node, Location, Host.rate, height_median, posterior)

cat("Mapping parent states...\n")
transitions <- df %>%
  filter(parent != node) %>%
  rename(child = node, child_Location = Location, child_Host = Host.rate, child_height = height_median, child_posterior = posterior) %>%
  left_join(lookup, by = c("parent" = "node")) %>%
  rename(parent_Location = Location, parent_Host = Host.rate)

cat("Filtering to actual transitions...\n")
jumps <- transitions %>%
  filter((parent_Location != child_Location) | (parent_Host != child_Host)) %>%
  mutate(
    from_state = paste(parent_Location, parent_Host, sep = " | "),
    to_state = paste(child_Location, child_Host, sep = " | ")
  )

cat("Summarizing transitions...\n")
summary_table <- jumps %>%
  group_by(parent_Location, parent_Host, child_Location, child_Host, from_state, to_state) %>%
  summarise(
    count = n(), 
    jump_tmrcas = paste(round(child_height, 2), collapse=", "),
    jump_posteriors = paste(round(child_posterior, 2), collapse=", "),
    .groups = "drop"
  )

# Append BF and PP
summary_table <- summary_table %>%
  rowwise() %>%
  mutate(
    Location_PP = if_else(parent_Location != child_Location, get_indicator_mean("Location", parent_Location, child_Location), NA_real_),
    Location_BF = if_else(!is.na(Location_PP), (Location_PP / (1 - Location_PP)) / prior_odds_loc, NA_real_),
    Host_PP = if_else(parent_Host != child_Host, get_indicator_mean("Host", parent_Host, child_Host), NA_real_),
    Host_BF = if_else(!is.na(Host_PP), (Host_PP / (1 - Host_PP)) / prior_odds_host, NA_real_)
  ) %>%
  ungroup() %>%
  # Handle Inf BFs
  mutate(
    Location_BF = ifelse(Location_PP == 1, Inf, Location_BF),
    Host_BF = ifelse(Host_PP == 1, Inf, Host_BF)
  )

# Filter to only Ecuador-related jumps and format output
summary_table <- summary_table %>%
  filter(grepl("ecuador", from_state, ignore.case=TRUE) | grepl("ecuador", to_state, ignore.case=TRUE)) %>%
  arrange(desc(count)) %>%
  select(from_state, to_state, count, jump_tmrcas, jump_posteriors, Location_PP, Location_BF, Host_PP, Host_BF)

write.csv(summary_table, out_file, row.names = FALSE)
cat(sprintf("Saved summary to %s\n", out_file))

