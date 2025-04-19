# Main Script
# This file runs the workflow for network meta-analysis and meta-regression

# Check environment setup
if (!file.exists("setup.R")) {
  stop("setup.R file not found. Please check if you are in the correct directory.")
}
source("setup.R")

# Load libraries and project sources
message("Loading libraries and project sources...")
library(here)  # for path management
library(yaml)  # for config loading
library(logger) # for structured logging
library(fs)    # for file operations

# Load functions
source(here("R", "01_utils.R"))
source(here("R", "02_data_prep.R"))
source(here("R", "03_models.R"))
source(here("R", "04_visualization.R"))

# Load configuration
config <- yaml::read_yaml(here("config", "config.yml"))
message("Configuration loaded successfully")

# Create timestamp directory for results
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_dir <- path(here("results"), timestamp)
dir_create(results_dir, recurse = TRUE)
message(paste("Results directory created:", results_dir))

# Update link to latest results
latest_link <- path(here("results"), "latest")
if (dir_exists(latest_link)) dir_delete(latest_link)
if (Sys.info()["sysname"] == "Windows") {
  # On Windows, create a batch file instead of symbolic link
  batch_content <- paste0("@echo off\ncd ", normalizePath(results_dir))
  writeLines(batch_content, path(here("results"), "latest.bat"))
  message("On Windows, created 'latest.bat' file.")
} else {
  # On Unix systems, create a symbolic link
  link_create(results_dir, latest_link)
  message("Symbolic link to latest results updated.")
}

# Logger setup
log_appender(appender_file(path(results_dir, "run.log")))
log_threshold(INFO)
log_info("Analysis started")

# Copy configuration file (for record)
file_copy(here("config", "config.yml"), path(results_dir, "config_used.yml"))

# Set random seed (for reproducibility)
set.seed(config$seed)
log_info(paste("Random seed set:", config$seed))

# Data preparation
log_info("Data preparation started")
# Use original data file with relative path
data_file_path <- "Sepsis_nma_data.csv"  # カレントディレクトリの直下に存在することを確認済み
message(paste("Using original data file:", data_file_path))
data <- prepare_data(
  file_path = data_file_path,
  markers = config$data$markers
)
saveRDS(data, path(results_dir, "prepared_data.rds"))
log_info("Data preparation completed")

# Run NMA model
log_info("Network meta-analysis started")
nma_results <- run_nma_model(
  data = data,
  iter = config$nma$iter,
  chains = config$nma$chains,
  cores = config$nma$cores,
  adapt_delta = config$nma$adapt_delta,
  stepsize = config$nma$stepsize,
  max_treedepth = config$nma$max_treedepth,
  results_dir = results_dir
)
saveRDS(nma_results, path(results_dir, "nma_results.rds"))

# Output NMA results to CSV
write.csv(nma_results$mu, path(results_dir, "nma_sensitivity_specificity.csv"))
write.csv(nma_results$dor, path(results_dir, "nma_dor.csv"))
write.csv(nma_results$sindex, path(results_dir, "nma_sindex.csv"))
log_info("Network meta-analysis completed")

# Run meta-regression model (if enabled in configuration)
regression_results <- NULL
if (config$meta_regression$run) {
  log_info("Meta-regression model started")
  regression_results <- run_regression_model(
    data = data,
    covariate = config$meta_regression$covariate,
    iter = config$meta_regression$iter,
    chains = config$meta_regression$chains,
    cores = config$meta_regression$cores,
    adapt_delta = config$meta_regression$adapt_delta,
    stepsize = config$meta_regression$stepsize,
    max_treedepth = config$meta_regression$max_treedepth,
    results_dir = results_dir
  )
  saveRDS(regression_results, path(results_dir, "regression_results.rds"))
  
  # Output meta-regression results to CSV
  write.csv(regression_results$mu, path(results_dir, "regression_sensitivity_specificity.csv"))
  write.csv(regression_results$dor, path(results_dir, "regression_dor.csv"))
  write.csv(regression_results$sindex, path(results_dir, "regression_sindex.csv"))
  write.csv(regression_results$beta1, path(results_dir, "regression_coefficients.csv"))
  log_info("Meta-regression model completed")
}

# Results visualization
log_info("Visualization creation started")
plots <- create_visualizations(
  nma_results = nma_results, 
  regression_results = regression_results, 
  reference_marker = config$data$reference_marker
)

# Save plots
for (i in seq_along(plots)) {
  plot_name <- names(plots)[i]
  plot_path <- path(results_dir, paste0(plot_name, ".png"))
  ggsave(plot_path, plots[[i]], width = 10, height = 8, dpi = 300)
  log_info(paste("Plot saved:", plot_name))
}
log_info("Visualization creation completed")

# Record execution information
run_info <- list(
  timestamp = Sys.time(),
  r_version = R.version.string,
  platform = Sys.info()["sysname"],
  packages = sapply(required_packages, function(pkg) packageVersion(pkg)),
  config = config,
  execution_time = list(
    nma = nma_results$execution_info$elapsed_time
  )
)

if (!is.null(regression_results)) {
  run_info$execution_time$regression <- regression_results$execution_info$elapsed_time
}

# Save as JSON and RDS
jsonlite::write_json(run_info, path(results_dir, "run_info.json"), pretty = TRUE, auto_unbox = TRUE)
saveRDS(run_info, path(results_dir, "run_info.rds"))
log_info("Execution information recorded")

# Delete old results (if enabled in configuration)
if (config$results$prune_old_results) {
  log_info("Deleting old results...")
  pruned_count <- prune_old_results(
    results_path = here("results"), 
    keep_last = config$results$keep_last_n_results
  )
  log_info(paste("Number of old result directories deleted:", pruned_count))
}

log_info("Analysis completed")
message("Analysis completed. Results saved to directory:")
message(results_dir)
