source("1_Helper_data_gen.R")
source("2_Helper_data_process.R")
source("3_Helper_run.R")

run_simulation <- function(n_replications = 25, run_number = 1) {
  library(pROC)
  
  # Create log file for diagnostics
  output_dir <- "simulation_results7"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  log_file <- file.path(output_dir, sprintf("diagnostics_run%d.txt", run_number))
  
  # Function to log messages
  log_msg <- function(...) {
    msg <- paste0(...)
    cat(msg, "\n")
    cat(msg, "\n", file = log_file, append = TRUE)
  }
  
  log_msg(paste(rep("=", 80), collapse = ""))
  log_msg("SIMULATION RUN ", run_number, " STARTED: ", Sys.time())
  log_msg("STRATEGY: BOTH (na.rm and epsilon)")
  log_msg(paste(rep("=", 80), collapse = ""))
  
  data_types <- c("LDA")
  scenarios <- c("nonrare_noninfo")
  
  log_msg("Starting simulation with ", n_replications, " replications")
  results_list_narm <- vector("list", n_replications)
  results_list_epsilon <- vector("list", n_replications)
  
  for (rep in 1:n_replications) {
    log_msg("\n", paste(rep("=", 60), collapse = ""))
    log_msg(sprintf("REPLICATION %d/%d", rep, n_replications))
    log_msg(paste(rep("=", 60), collapse = ""))
    
    # Initialize results dataframes for BOTH strategies
    rep_results_narm <- data.frame(
      replication = integer(),
      data_type = character(),
      scenario = character(),
      algorithm = character(),
      set = character(),
      n_total = integer(),
      n_included = integer(),
      n_missing = integer(),
      n_na_predictions = integer(),
      n_na_true_values = integer(),
      n_incomplete_cases = integer(),
      mse = numeric(),
      auc = numeric(),
      f1 = numeric(),
      precision = numeric(),
      recall = numeric(),
      accuracy = numeric(),
      prob_mse = numeric(),
      prob_mae = numeric(),
      prob_correlation = numeric(),
      optimal_threshold = numeric(),
      stringsAsFactors = FALSE
    )
    
    rep_results_epsilon <- rep_results_narm  # Same structure
    
    for (data_type in data_types) {
      for (scenario in scenarios) {
        log_msg(sprintf("\n  Processing %s - %s", data_type, scenario))
        
        # Generate and process data ONCE
        log_msg("    Generating data...")
        dataset <- switch(data_type,
                          "ideal" = create_simulation_data(scenario),
                          "ori" = create_complex_data(n = 10000, scenario = scenario),
                          "LDA" = create_LDA_data(n = 10000, scenario = scenario))
        
        log_msg("    Processing data...")
        processed_phenorm <- process_phenorm(dataset, data_type, scenario, valid_label = "test_set")
        processed_map <- process_map(dataset, data_type, scenario)
        processed_surelda <- process_surelda(dataset, data_type, scenario)
        processed_surelda_map <- process_surelda_map(dataset, data_type, scenario)
        log_msg("    done")
        
        # Run algorithms ONCE
        log_msg("    Running algorithms...")
        phenorm_results <- run_phenorm_model(processed_phenorm)
        map_results <- run_map_model(processed_map)
        surelda_results <- run_surelda_model(processed_surelda)
        surelda_phenorm_results <- run_surelda_phenorm(processed_surelda, phenorm_results)
        surelda_map_results <- run_surelda_map(processed_surelda_map)
        surelda_map_p_results <- run_surelda_map_p(processed_surelda_map, map_results)
        icd_logit_results <- run_icd_logit(dataset)
        log_msg("    done")
        
        log_msg("\n    Calculating metrics with BOTH strategies...")
        
        # Helper function to process results for both strategies
        process_algorithm_results <- function(algo_name, results_obj, 
                                              pred_train, pred_test,
                                              idx_train, idx_test) {
          if (!is.null(results_obj)) {
            for (set_type in c("train", "test")) {
              predictions <- if(set_type == "train") pred_train else pred_test
              indices <- if(set_type == "train") idx_train else idx_test
              
              if (!is.null(predictions) && !is.null(indices)) {
                true_values <- dataset$probability[indices]
                true_outcomes <- dataset$outcome[indices]
                
                # Calculate metrics with BOTH strategies
                metrics_narm <- calculate_metrics(
                  predictions, true_values, true_outcomes,
                  n_total = length(indices),
                  n_included = length(predictions),
                  use_epsilon = FALSE  # na.rm strategy
                )
                
                metrics_epsilon <- calculate_metrics(
                  predictions, true_values, true_outcomes,
                  n_total = length(indices),
                  n_included = length(predictions),
                  use_epsilon = TRUE  # epsilon strategy
                )
                
                log_msg(sprintf("      %s (%s) - narm: MSE=%.4f, AUC=%.4f, incomplete=%d",
                                algo_name, set_type, metrics_narm$mse, metrics_narm$auc, 
                                metrics_narm$n_incomplete_cases))
                log_msg(sprintf("      %s (%s) - epsilon: MSE=%.4f, AUC=%.4f, incomplete=%d",
                                algo_name, set_type, metrics_epsilon$mse, metrics_epsilon$auc, 
                                metrics_epsilon$n_incomplete_cases))
                
                # Store results for both strategies
                rep_results_narm <<- rbind(rep_results_narm, data.frame(
                  replication = rep,
                  data_type = data_type,
                  scenario = scenario,
                  algorithm = algo_name,
                  set = set_type,
                  metrics_narm,
                  stringsAsFactors = FALSE
                ))
                
                rep_results_epsilon <<- rbind(rep_results_epsilon, data.frame(
                  replication = rep,
                  data_type = data_type,
                  scenario = scenario,
                  algorithm = algo_name,
                  set = set_type,
                  metrics_epsilon,
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
        }
        
        # Process all algorithms
        process_algorithm_results("phenorm", phenorm_results,
                                  phenorm_results$train_preds$Aggregate,
                                  phenorm_results$test_preds$Aggregate,
                                  phenorm_results$train_idx,
                                  phenorm_results$test_idx)
        
        process_algorithm_results("map", map_results,
                                  map_results$train$scores,
                                  map_results$test$scores,
                                  map_results$train_idx,
                                  map_results$test_idx)
        
        process_algorithm_results("surelda", surelda_results,
                                  surelda_results$train$probs,
                                  surelda_results$test$probs,
                                  surelda_results$train_idx,
                                  surelda_results$test_idx)
        
        process_algorithm_results("surelda_phenorm", surelda_phenorm_results,
                                  surelda_phenorm_results$train$probs,
                                  surelda_phenorm_results$test$probs,
                                  surelda_phenorm_results$train_idx,
                                  surelda_phenorm_results$test_idx)
        
        process_algorithm_results("surelda_map", surelda_map_results,
                                  surelda_map_results$train$probs,
                                  surelda_map_results$test$probs,
                                  surelda_map_results$train_idx,
                                  surelda_map_results$test_idx)
        
        # Process priors
        if (!is.null(surelda_results) &&
            !is.null(surelda_results$train$prior) &&
            !is.null(surelda_results$test$prior)) {
          process_algorithm_results("surelda_prior", surelda_results,
                                    surelda_results$train$prior,
                                    surelda_results$test$prior,
                                    surelda_results$train_idx,
                                    surelda_results$test_idx)
        }
        
        if (!is.null(surelda_map_results) &&
            !is.null(surelda_map_results$train$prior) &&
            !is.null(surelda_map_results$test$prior)) {
          process_algorithm_results("surelda_map_prior", surelda_map_results,
                                    surelda_map_results$train$prior,
                                    surelda_map_results$test$prior,
                                    surelda_map_results$train_idx,
                                    surelda_map_results$test_idx)
        }
        
        process_algorithm_results("surelda_map_p", surelda_map_p_results,
                                  surelda_map_p_results$train$probs,
                                  surelda_map_p_results$test$probs,
                                  surelda_map_p_results$train_idx,
                                  surelda_map_p_results$test_idx)
        
        process_algorithm_results("icd_logit", icd_logit_results,
                                  icd_logit_results$train$probs,
                                  icd_logit_results$test$probs,
                                  icd_logit_results$train_idx,
                                  icd_logit_results$test_idx)
        
        log_msg("    done")
      }
    }
    
    results_list_narm[[rep]] <- rep_results_narm
    results_list_epsilon[[rep]] <- rep_results_epsilon
  }
  
  # Combine results for both strategies
  log_msg("\n", paste(rep("=", 60), collapse = ""))
  log_msg("COMBINING RESULTS")
  log_msg(paste(rep("=", 60), collapse = ""))
  
  results_narm <- do.call(rbind, results_list_narm)
  results_epsilon <- do.call(rbind, results_list_epsilon)
  
  # Calculate summary statistics for both
  summary_stats_narm <- results_narm %>%
    group_by(data_type, scenario, algorithm, set) %>%
    summarise(across(c(n_total, n_included, n_missing, 
                       n_na_predictions, n_na_true_values, n_incomplete_cases,
                       mse, auc, f1, precision, recall, accuracy, 
                       prob_mse, prob_mae, prob_correlation, optimal_threshold),
                     list(mean = mean, sd = sd)),
              .groups = 'drop')
  
  summary_stats_epsilon <- results_epsilon %>%
    group_by(data_type, scenario, algorithm, set) %>%
    summarise(across(c(n_total, n_included, n_missing, 
                       n_na_predictions, n_na_true_values, n_incomplete_cases,
                       mse, auc, f1, precision, recall, accuracy, 
                       prob_mse, prob_mae, prob_correlation, optimal_threshold),
                     list(mean = mean, sd = sd)),
              .groups = 'drop')
  
  log_msg("\nSaving results...")
  save(results_narm, summary_stats_narm,
       file = file.path(output_dir,
                        sprintf("simulation_results_narm_run%d.RData", run_number)))
  
  save(results_epsilon, summary_stats_epsilon,
       file = file.path(output_dir,
                        sprintf("simulation_results_epsilon_run%d.RData", run_number)))
  
  log_msg("\n", paste(rep("=", 80), collapse = ""))
  log_msg("SIMULATION COMPLETED: ", Sys.time())
  log_msg("Results saved to: ", output_dir)
  log_msg(paste(rep("=", 80), collapse = ""))
  
  return(list(
    narm = list(detailed = results_narm, summary = summary_stats_narm),
    epsilon = list(detailed = results_epsilon, summary = summary_stats_epsilon)
  ))
}

# Calculate metrics function WITH strategy parameter
calculate_metrics <- function(predictions, true_values, true_outcomes,
                              n_total, n_included, 
                              use_epsilon = TRUE) {
  n_missing <- n_total - n_included
  epsilon <- 1e-10
  
  # Track NAs BEFORE any modifications
  n_na_predictions_original <- sum(is.na(predictions))
  n_na_true_values_original <- sum(is.na(true_values))
  n_incomplete_cases_original <- sum(is.na(predictions) | is.na(true_values))
  
  # Apply strategy
  if(use_epsilon) {
    # EPSILON: Replace NAs and ensure minimum value
    true_values <- ifelse(is.na(true_values), epsilon, pmax(true_values, epsilon))
  }
  
  predictions <- as.numeric(predictions)
  if(any(is.na(predictions))) {
    if (!use_epsilon) {
      # NA.RM: Impute with mean
      predictions[is.na(predictions)] <- mean(predictions, na.rm = TRUE)
    } else {
      # EPSILON: Replace with epsilon
      predictions[is.na(predictions)] <- epsilon
    }
  }
  
  true_outcomes <- as.numeric(true_outcomes)
  
  auc <- precision <- recall <- f1 <- accuracy <- NA
  prob_mse <- prob_mae <- prob_correlation <- optimal_threshold <- NA
  
  tryCatch({
    roc_obj <- roc(true_outcomes, predictions, quiet = TRUE)
    auc <- as.numeric(auc(roc_obj))
    
    coords_result <- coords(roc_obj, "best", best.method = "youden", quiet = TRUE)
    optimal_threshold <- coords_result$threshold
    
    binary_preds <- predictions >= optimal_threshold
    
    cm <- table(Predicted = binary_preds, Actual = true_outcomes)
    
    if (nrow(cm) == 2 && ncol(cm) == 2) {
      precision <- cm[2,2] / sum(cm[2,])
      recall <- cm[2,2] / sum(cm[,2])
      
      if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
        f1 <- 2 * (precision * recall) / (precision + recall)
      }
      
      accuracy <- sum(diag(cm)) / sum(cm)
    }
    
    # Set na.rm based on strategy
    use_narm <- !use_epsilon
    prob_mse <- mean((predictions - true_values)^2, na.rm = use_narm)
    prob_mae <- mean(abs(predictions - true_values), na.rm = use_narm)
    prob_correlation <- cor(predictions, true_values, use = "complete.obs")
    
  }, error = function(e) {
    warning("Error in metrics calculation: ", e$message)
  })
  
  data.frame(
    n_total = n_total,
    n_included = n_included,
    n_missing = n_missing,
    n_na_predictions = n_na_predictions_original,
    n_na_true_values = n_na_true_values_original,
    n_incomplete_cases = n_incomplete_cases_original,
    mse = prob_mse,
    auc = auc,
    f1 = f1,
    precision = precision,
    recall = recall,
    accuracy = accuracy,
    prob_mse = prob_mse,
    prob_mae = prob_mae,
    prob_correlation = prob_correlation,
    optimal_threshold = optimal_threshold
  )
}

library(pROC)
args <- commandArgs(trailingOnly = TRUE)
run_number <- as.numeric(args[1])
if (is.na(run_number)) {
  stop("Run number must be provided as a command line argument")
}
set.seed(20250606 + run_number * 1000)
results <- run_simulation(n_replications = 32, run_number = run_number)