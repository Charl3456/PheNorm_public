library(dplyr)
library(Matrix)
library(logger)
library(PheNorm)
library(MAP)
library(sureLDA)

# Data split----------------------------------------------------------------------------------------------------------
split_data <- function(data, test_size = 200, seed = 1234) {
  if (!is.numeric(test_size) || test_size < 1) {
    stop("test_size must be a positive number")
  }

  total_rows <- nrow(data)
  if (test_size >= total_rows) {
    stop("test_size must be less than total number of rows")
  }

  # set.seed(seed)
  test_indices <- sample(1:total_rows, test_size)
  train_indices <- setdiff(1:total_rows, test_indices)

  list(
    train = data[train_indices,],
    test = data[test_indices,],
    train_idx = train_indices,
    test_idx = test_indices
  )
}

# Helper functions ---------------------------------------------------------------------------------------------------
predict.PheNorm <- function(phenorm_model = NULL, newdata = NULL,
                            silver_labels = NULL, features = NULL,
                            utilization = NULL, use_empirical_sd = TRUE,
                            aggregate_labels = silver_labels,
                            na.rm = TRUE) {
  # Input validation
  if (is.null(phenorm_model)) stop("Please enter a fitted PheNorm model.")
  if (is.null(newdata)) stop("Please enter a new dataset for predictions.")
  if (is.null(silver_labels)) stop("silver_labels must be specified")

  # Ensure newdata is a matrix
  newmat <- as.matrix(newdata)

  # Validate column existence
  missing_cols <- setdiff(silver_labels, colnames(newmat))
  if (length(missing_cols) > 0) {
    stop("Missing columns in newdata: ", paste(missing_cols, collapse = ", "))
  }

  # Extract and normalize silver labels
  silver_label_data <- newmat[, silver_labels, drop = FALSE]
  util <- newmat[, utilization]

  # Handle potential NA values in util
  if (any(is.na(util))) {
    warning("NA values found in utilization, replacing with mean")
    util[is.na(util)] <- mean(util, na.rm = TRUE)
  }

  # Normalize silver labels
  normalized_silver_labels <- silver_label_data -
    PheNorm:::VTM(phenorm_model$alpha, nrow(newdata)) * util

  # Process features if present
  if (!is.null(features)) {
    feature_data <- newmat[, features, drop = FALSE]
    normalized_data <- cbind(normalized_silver_labels, feature_data)
  } else {
    normalized_data <- normalized_silver_labels
  }

  # Calculate phenorm scores
  phenorm_score <- try(normalized_data %*% as.matrix(phenorm_model$betas))
  if (inherits(phenorm_score, "try-error")) {
    stop("Error in computing phenorm scores: ", phenorm_score)
  }

  # Process each column separately for better error handling
  posterior_probs <- matrix(NA, nrow = nrow(phenorm_score), ncol = ncol(phenorm_score))
  colnames(posterior_probs) <- colnames(phenorm_score)

  for (i in 1:ncol(phenorm_score)) {
    tryCatch({
      fit <- PheNorm:::normalmixEM2comp2(
        phenorm_score[, i],
        lambda = 0.5,
        mu = quantile(phenorm_model$scores[, i], probs = c(1/3, 2/3), na.rm = na.rm),
        sigsqrd = ifelse(use_empirical_sd,
                         max(sd(phenorm_model$scores[, i], na.rm = na.rm) / 2, 0.1),
                         1)
      )
      posterior_probs[, i] <- fit$posterior[, 2]
    }, error = function(e) {
      warning(sprintf("Error processing column %d: %s", i, e$message))
      posterior_probs[, i] <- NA
    })
  }

  # Calculate aggregate probabilities
  posterior_probs_vote <- rowMeans(posterior_probs[, aggregate_labels, drop = FALSE],
                                   na.rm = TRUE)

  all_posterior_probs <- cbind.data.frame(
    posterior_probs,
    "Aggregate" = posterior_probs_vote
  )

  return(all_posterior_probs)
}

## Run Phenorm Helper functions ---------------------------------------------------------------------------------------------------------------------------------------------------
run_phenorm <- function(train = NULL, test = NULL, silver_labels = "", features = "",
                        utilization = "", weight = "", seed = 1234,
                        aggregate_labels = silver_labels, ...) {
  # Input validation
  if (is.null(train) || is.null(test)) {
    stop("Both train and test datasets must be provided")
  }

  # Get additional parameters
  params <- list(...)
  corrupt_rate <- params$corrupt.rate %||% 0.3
  train_size <- params$train.size %||% (10 * nrow(train))
  cat("Running PheNorm model...\n")

  # Fit model with error handling
  phenorm_fit <- tryCatch({
    # set.seed(seed)
    phenorm_prob(
      nm.logS.ori = silver_labels,
      nm.utl = utilization,
      nm.wt = weight,
      dat = train,
      nm.X = features,
      corrupt.rate = corrupt_rate,
      train.size = train_size
    )
  }, error = function(e) {
    stop("Error in PheNorm fitting: ", e$message)
    return(NULL)
  })

  if (is.null(phenorm_fit)) {
    return(NULL)
  }

  # Generate predictions with error handling
  train_preds <- tryCatch({
    # set.seed(seed)
    predict.PheNorm(
      phenorm_model = phenorm_fit,
      newdata = train,
      silver_labels = silver_labels,
      features = features,
      utilization = utilization,
      aggregate_labels = aggregate_labels
    )
  }, error = function(e) {
    warning("Error in training predictions: ", e$message)
    return(NULL)
  })

  test_preds <- tryCatch({
    predict.PheNorm(
      phenorm_model = phenorm_fit,
      newdata = test,
      silver_labels = silver_labels,
      features = features,
      utilization = utilization,
      aggregate_labels = aggregate_labels
    )
  }, error = function(e) {
    warning("Error in test predictions: ", e$message)
    return(NULL)
  })
cat("PheNorm done \n")
  return(list(
    "fit" = phenorm_fit,
    "train_preds" = train_preds,
    "test_preds" = test_preds
  ))
}

# Add create_summary_tables function
create_summary_tables <- function(summary_stats) {
  # Group by algorithm and set
  tables <- list()

  # Split summary stats by metric type
  metrics <- c("mse", "auc", "f1", "precision", "recall", "accuracy",
               "prob_mse", "prob_mae", "prob_correlation")

  for (metric in metrics) {
    mean_col <- paste0(metric, "_mean")
    sd_col <- paste0(metric, "_sd")

    # Create table for current metric
    metric_table <- summary_stats %>%
      select(data_type, scenario, algorithm, set, !!mean_col, !!sd_col) %>%
      mutate(
        value = sprintf("%.3f (%.3f)", get(mean_col), get(sd_col))
      ) %>%
      select(data_type, scenario, algorithm, set, value) %>%
      pivot_wider(
        names_from = c(set),
        values_from = value
      )

    tables[[metric]] <- metric_table
  }

  return(tables)
}

  phenorm_prob <- function(nm.logS.ori, nm.utl, nm.wt, dat, nm.X = NULL, corrupt.rate = 0.3, train.size = 10 * nrow(dat)) {
    dat <- as.matrix(dat)
    S.ori <- dat[, nm.logS.ori, drop = FALSE]
    utl <- dat[, nm.utl, drop = FALSE]
    if (nm.wt == "") {
      wt <- NULL
    } else {
      wt <- dat[, nm.wt, drop = FALSE]
    }
    a.hat <- apply(S.ori, 2, function(S) {PheNorm:::findMagicNumber(S, utl)$coef})
    S.norm <- S.ori - PheNorm:::VTM(a.hat, nrow(dat)) * as.vector(utl)
    if (!is.null(nm.X)) {
      X <- as.matrix(dat[, nm.X, drop = FALSE])
      if (length(unique(utl)) == 1) {
        SX.norm <- cbind(S.norm, X, wt)
      } else {
        SX.norm <- cbind(S.norm, X, utl, wt)
      }
      id <- sample(1:nrow(dat), train.size, replace = TRUE)
      SX.norm.corrupt <- apply(SX.norm[id, ], 2,
                               function(x) {ifelse(rbinom(length(id), 1, corrupt.rate), mean(x), x)}
      )
      sx_norm_df <- as.data.frame(SX.norm.corrupt)
      if (nm.wt == "") {
        weights <- rep(1, nrow(sx_norm_df))
      } else {
        weights <- sx_norm_df[[nm.wt]]
      }
      # seperate here
      b.all <- apply(S.norm, 2, function(ss) {
        lm(ss[id] ~ . - 1, data = sx_norm_df[, !(names(sx_norm_df) %in% nm.wt)],
           weights = weights)$coef
      })
      b.all[is.na(b.all)] <- 0
      S.norm <- as.matrix(SX.norm[, !(colnames(SX.norm) %in% nm.wt)]) %*% b.all
      if (length(unique(utl)) > 1) {
        b.all <- b.all[-dim(b.all)[1], ]
      }
    } else {
      b.all <- NULL
    }
    if (length(nm.logS.ori) > 1) {
      postprob <- apply(S.norm, 2,
                        function(x) {
                          fit = PheNorm:::normalmixEM2comp2(x, lambda = 0.5,
                                                            mu = quantile(x, probs=c(1/3, 2/3)), sigsqrd = 1
                          )
                          fit$posterior[, 2]
                        }
      )
      list("probs" = rowMeans(postprob, na.rm = TRUE), "prob_individual" = postprob, "betas" = b.all, "scores" = S.norm,
           "alpha" = a.hat)

    } else {
      fit <- PheNorm:::normalmixEM2comp2(unlist(S.norm), lambda = 0.5,
                                         mu = quantile(S.norm, probs=c(1/3, 2/3)), sigsqrd = 1
      )
      list("probs" = fit$posterior[,2], "prob_individual" = postprob, "betas" = b.all, "scores" = S.norm,
           "alpha" = a.hat)
    }
  }

predict.PheNorm <- function(phenorm_model = NULL, newdata = NULL,
                            silver_labels = NULL, features = NULL,
                            utilization = NULL, use_empirical_sd = TRUE,
                            aggregate_labels = silver_labels,
                            na.rm = TRUE) {
  if (is.null(phenorm_model)) {
    stop("Please enter a fitted PheNorm model.")
  }
  if (is.null(newdata)) {
    stop("Please enter a new dataset to get PheNorm predictions on.")
  }
  # normal mixture normalization on the new dataset
  newmat <- as.matrix(newdata)
  silver_label_data <- newmat[, silver_labels, drop = FALSE]
  util <- newmat[, utilization]
  normalized_silver_labels <- silver_label_data -
    PheNorm:::VTM(phenorm_model$alpha, nrow(newdata)) * util
  # also add features if we used them
  if (!is.null(features)) {
    feature_data <- newmat[, features, drop = FALSE]
    normalized_data <- cbind(normalized_silver_labels, feature_data)
  } else {
    normalized_data <- normalized_silver_labels
  }
  # predict, for each silver label separately and overall ("voting")
  phenorm_score <- normalized_data %*% as.matrix(phenorm_model$betas)
  original_phenorm_score <- phenorm_model$scores
  posterior_probs <- do.call(cbind, sapply(1:ncol(phenorm_score), function(i) {
    fit <- PheNorm:::normalmixEM2comp2(phenorm_score[, i], lambda = 0.5,
                                       mu = quantile(original_phenorm_score[, i], probs = c(1/3, 2/3), na.rm = na.rm),
                                       sigsqrd = ifelse(use_empirical_sd, sd(original_phenorm_score[, i], na.rm = na.rm) / 2, 1))
    fit$posterior[, 2]
  }, simplify = FALSE))
  colnames(posterior_probs) <- colnames(phenorm_score)
  posterior_probs_vote <- rowMeans(posterior_probs[, aggregate_labels])
  all_posterior_probs <- cbind.data.frame(posterior_probs, "Aggregate" = posterior_probs_vote)
  return(all_posterior_probs)
}


##Run Phenorm----------------------------------------------------------------------------------------------------------------------------
run_phenorm_model <- function(analysis_data,
                              utilization_var="Utiliz",
                              weight_var="Sampling_Weight",
                              corrupt_rate=0.3,
                              train_size_mult=13,
                              seed=4747) {

  # Split data into train and test
  # split_result <- split_data(analysis_data$test, test_size = 200, seed = seed)
  train_data <- analysis_data$train
  test_data <- analysis_data$test

  silver_labels <- analysis_data$silver_labels
  features <- setdiff(names(train_data),
                      c(silver_labels, utilization_var, weight_var))

  tryCatch({
    phenorm_analysis <- run_phenorm(
      train = train_data,
      test = test_data,
      silver_labels = silver_labels,
      aggregate_labels = silver_labels,
      features = features,
      utilization = utilization_var,
      weight = weight_var,
      corrupt.rate = corrupt_rate,
      train.size = train_size_mult * nrow(train_data),
      seed = seed
    )
    phenorm_analysis$train_idx <- which(analysis_data$all$test_set == 0)  # Add these lines
    phenorm_analysis$test_idx <- which(analysis_data$all$test_set == 1)
    return(phenorm_analysis)

  }, error = function(e) {
    warning(paste("Error in PheNorm analysis:", e$message))
    return(NULL)
  })
}

## Run MAP - Updated to extract train/test split
run_map_model <- function(data, params = list(
  verbose = TRUE,
  yes.con = FALSE,
  full.output = TRUE
)) {
  # Input validation
  if (is.null(data$mat) || is.null(data$note)) {
    stop("data must contain both 'mat' and 'note' components")
  }

  # Extract train/test split from processed data
  train_idx <- which(data$test_set == 0)
  test_idx <- which(data$test_set == 1)

  # Prepare split data and ensure sparse matrix format
  train_mat <- as(data$mat[train_idx,], "sparseMatrix")
  train_note <- Matrix(as.matrix(data$note[train_idx,]), sparse=TRUE)
  test_mat <- as(data$mat[test_idx,], "sparseMatrix")
  test_note <- Matrix(as.matrix(data$note[test_idx,]), sparse=TRUE)

  # Define custom fitproj_flexmix that uses MAP's fit_flexmix
  custom_fitproj_flexmix <- function(tmpfm, note, family, tmpfm2, dat.tmp, sample_size) {
    pkg_fit_flexmix <- getFromNamespace("fit_flexmix", "MAP")

    # Keep original split
    train_samples <- seq_len(nrow(dat.tmp))

    m_data = dat.tmp
    m_note = note
    rownames(m_data) = train_samples

    # Use all training samples
    m_fit = m_data
    m_fit_note = m_note

    fit_obj = pkg_fit_flexmix(tmpfm, m_fit_note, family, tmpfm2, m_fit)

    # Get cluster values using methods instead of direct slot access
    if(methods::is(fit_obj$tmpfit, "flexmix")) {
      cluster_values <- flexmix::clusters(fit_obj$tmpfit)
    } else {
      cluster_values <- fit_obj$tmpfit$cluster
    }

    # Return training data and results
    list(m_data = m_data,
         m_post = flexmix::posterior(fit_obj$tmpfit),
         cluster = cluster_values)
  }

  # Get original fitproj_flexmix if it exists
  orig_fitproj_flexmix <- if(exists("fitproj_flexmix", envir = asNamespace("MAP"))) {
    get("fitproj_flexmix", envir = asNamespace("MAP"))
  } else {
    NULL
  }

  # Replace fitproj_flexmix in MAP namespace
  if(!environmentIsLocked(asNamespace("MAP"))) {
    unlockBinding("fitproj_flexmix", asNamespace("MAP"))
    assign("fitproj_flexmix", custom_fitproj_flexmix, asNamespace("MAP"))
    lockBinding("fitproj_flexmix", asNamespace("MAP"))
  }

  # Run MAP on training data
  cat("Running MAP on training data...\n")
  train_results <- tryCatch({
    MAP::MAP(
      mat = train_mat,
      note = train_note,
      yes.con = params$yes.con,
      full.output = TRUE,
      subset_sample = FALSE,  # Use all training data
      verbose = params$verbose
    )
  }, error = function(e) {
    cat("Error in MAP training:", e$message, "\n")
    if(!is.null(e$call)) cat("Call:", deparse(e$call), "\n")
    return(NULL)
  })

  # Run MAP on test data
  cat("Running MAP on test data...\n")
  test_results <- tryCatch({
    MAP::MAP(
      mat = test_mat,
      note = test_note,
      yes.con = params$yes.con,
      full.output = TRUE,
      subset_sample = FALSE,  # Use all test data
      verbose = params$verbose
    )
  }, error = function(e) {
    cat("Error in MAP testing:", e$message, "\n")
    if(!is.null(e$call)) cat("Call:", deparse(e$call), "\n")
    return(NULL)
  }, finally = {
    # Restore original fitproj_flexmix if it existed
    if(!is.null(orig_fitproj_flexmix) && !environmentIsLocked(asNamespace("MAP"))) {
      unlockBinding("fitproj_flexmix", asNamespace("MAP"))
      assign("fitproj_flexmix", orig_fitproj_flexmix, asNamespace("MAP"))
      lockBinding("fitproj_flexmix", asNamespace("MAP"))
    }
  })

  if(is.null(train_results) || is.null(test_results)) {
    cat("MAP analysis failed\n")
    return(NULL)
  }

  list(
    train = list(scores = train_results$scores),
    test = list(scores = test_results$scores),
    train_idx = train_idx,
    test_idx = test_idx,
    cut.MAP = train_results$cut.MAP
  )
}

## Run SureLDA - Updated to extract train/test split
run_surelda_model <- function(data, params = list(
  prior = "PheNorm",
  weight = "beta",
  nEmpty = 20,
  alpha = 100,
  beta = 100,
  burnin = 50,
  ITER = 150,
  verbose = TRUE
)) {
  cat("Starting sureLDA analysis...\n")

  # Extract train/test split from processed data
  train_idx <- which(data$test_set == 0)
  test_idx <- which(data$test_set == 1)

  train_data <- list(
    X = as.matrix(data$X[train_idx,]),
    ICD = as.matrix(data$ICD[train_idx,]),
    NLP = as.matrix(data$NLP[train_idx,]),
    HU = as.vector(data$HU[train_idx]),
    filter = as.matrix(data$filter[train_idx,])
  )

  test_data <- list(
    X = as.matrix(data$X[test_idx,]),
    ICD = as.matrix(data$ICD[test_idx,]),
    NLP = as.matrix(data$NLP[test_idx,]),
    HU = as.vector(data$HU[test_idx]),
    filter = as.matrix(data$filter[test_idx,])
  )

  tryCatch({
    train_results <- sureLDA(
      X = train_data$X,
      ICD = train_data$ICD,
      NLP = train_data$NLP,
      HU = train_data$HU,
      filter = train_data$filter,
      prior = params$prior,
      weight = params$weight,
      nEmpty = params$nEmpty,
      alpha = params$alpha,
      beta = params$beta,
      burnin = params$burnin,
      ITER = params$ITER,
      verbose = params$verbose
    )

    # Test using trained phi
    test_results <- sureLDA(
      X = test_data$X,
      ICD = test_data$ICD,
      NLP = test_data$NLP,
      HU = test_data$HU,
      filter = test_data$filter,
      prior = params$prior,
      weight = params$weight,
      nEmpty = params$nEmpty,
      phi = train_results$phi,
      verbose = params$verbose
    )
    cat("sureLDA done \n")
    list(
      train = train_results,
      test = test_results,
      train_idx = train_idx,
      test_idx = test_idx
    )
  }, error = function(e) {
    warning(paste("Error in sureLDA analysis:", e$message))
    return(NULL)
  })
}

##Run sureLDA with external PheNorm prior-----------------------------------------------------------------------------------------
run_surelda_phenorm <- function(data, phenorm_prior, params = NULL) {
  if (is.null(params)) {
    params <- list(
      weight = "uniform",
      nEmpty = 20,
      alpha = 100,
      beta = 100,
      burnin = 50,
      ITER = 150,
      verbose = TRUE
    )
  }
  
  cat("Running sureLDA with PheNorm prior...\n")

  # Comprehensive input validation
  if (is.null(data$HU)) {
    warning("Healthcare Utilization (HU) data is missing")
    return(NULL)
  }
  
  if (is.null(phenorm_prior$train_preds) || is.null(phenorm_prior$test_preds)) {
    warning("PheNorm predictions missing")
    return(NULL)
  }
  
  if (is.null(phenorm_prior$train_preds$Aggregate) || is.null(phenorm_prior$test_preds$Aggregate)) {
    warning("PheNorm Aggregate predictions missing")
    return(NULL)
  }

  # Extract train/test indices consistently
  train_indices <- which(data$test_set == 0)
  test_indices <- which(data$test_set == 1)
  
  # Validate dimensions
  expected_train_size <- length(train_indices)
  expected_test_size <- length(test_indices)
  actual_train_preds <- length(phenorm_prior$train_preds$Aggregate)
  actual_test_preds <- length(phenorm_prior$test_preds$Aggregate)
  
  if (expected_train_size != actual_train_preds) {
    warning(sprintf("Train size mismatch: expected %d, got %d predictions", 
                    expected_train_size, actual_train_preds))
    return(NULL)
  }
  
  if (expected_test_size != actual_test_preds) {
    warning(sprintf("Test size mismatch: expected %d, got %d predictions", 
                    expected_test_size, actual_test_preds))
    return(NULL)
  }

  # Prepare data
  train_data <- list(
    X = as.matrix(data$X[train_indices,]),
    ICD = as.matrix(data$ICD[train_indices,]),
    NLP = as.matrix(data$NLP[train_indices,]),
    HU = as.numeric(data$HU[train_indices]),
    filter = as.matrix(data$filter[train_indices,])
  )

  test_data <- list(
    X = as.matrix(data$X[test_indices,]),
    ICD = as.matrix(data$ICD[test_indices,]),
    NLP = as.matrix(data$NLP[test_indices,]),
    HU = as.numeric(data$HU[test_indices]),
    filter = as.matrix(data$filter[test_indices,])
  )

  # Extract and validate priors
  train_prior <- as.numeric(phenorm_prior$train_preds$Aggregate)
  test_prior <- as.numeric(phenorm_prior$test_preds$Aggregate)
  
  # Check for invalid values
  if (any(is.na(train_prior)) || any(is.na(test_prior))) {
    warning("NA values found in PheNorm priors")
    return(NULL)
  }
  
  # Convert to matrices and ensure valid probability range
  train_prior <- as.matrix(pmax(pmin(train_prior, 0.999), 0.001))
  test_prior <- as.matrix(pmax(pmin(test_prior, 0.999), 0.001))
  
  cat(sprintf("Prior validation: Train [%.3f, %.3f], Test [%.3f, %.3f]\n",
              min(train_prior), max(train_prior), min(test_prior), max(test_prior)))

  # Run sureLDA with comprehensive error handling
  tryCatch({
    train_results <- sureLDA(
      X = train_data$X,
      ICD = train_data$ICD,
      NLP = train_data$NLP,
      HU = train_data$HU,
      filter = train_data$filter,
      prior = train_prior,
      weight = params$weight,
      nEmpty = params$nEmpty,
      alpha = params$alpha,
      beta = params$beta,
      burnin = params$burnin,
      ITER = params$ITER,
      verbose = params$verbose
    )

    test_results <- sureLDA(
      X = test_data$X,
      ICD = test_data$ICD,
      NLP = test_data$NLP,
      HU = test_data$HU,
      filter = test_data$filter,
      prior = test_prior,
      weight = params$weight,
      nEmpty = params$nEmpty,
      phi = train_results$phi,
      verbose = params$verbose
    )

    return(list(
      train = train_results,
      test = test_results,
      train_idx = train_indices,
      test_idx = test_indices
    ))
    
  }, error = function(e) {
    warning(paste("Error in sureLDA with PheNorm prior:", e$message))
    cat("Train data dimensions:\n")
    cat("X:", dim(train_data$X), "\n")
    cat("ICD:", dim(train_data$ICD), "\n") 
    cat("NLP:", dim(train_data$NLP), "\n")
    cat("HU length:", length(train_data$HU), "\n")
    cat("filter:", dim(train_data$filter), "\n")
    cat("prior:", dim(train_prior), "\n")
    return(NULL)
  })
}

## Run sureLDA with MAP prior - Updated to extract train/test split
run_surelda_map <- function(data, params = list(
  weight = "uniform",
  nEmpty = 20,
  alpha = 100,
  beta = 100,
  burnin = 50,
  ITER = 150,
  verbose = TRUE
)) {
  cat("Starting sureLDA with MAP prior...\n")

  # Extract train/test split from processed data
  train_idx <- which(data$test_set == 0)
  test_idx <- which(data$test_set == 1)

  # Prepare train/test data
  train_data <- list(
    X = as.matrix(data$X[train_idx,]),
    ICD = as.matrix(data$ICD[train_idx,]),
    NLP = as.matrix(data$NLP[train_idx,]),
    HU = as.vector(data$HU[train_idx]),
    filter = as.matrix(data$filter[train_idx,])
  )

  test_data <- list(
    X = as.matrix(data$X[test_idx,]),
    ICD = as.matrix(data$ICD[test_idx,]),
    NLP = as.matrix(data$NLP[test_idx,]),
    HU = as.vector(data$HU[test_idx]),
    filter = as.matrix(data$filter[test_idx,])
  )

  tryCatch({
    # Train with internal MAP prior
    train_results <- sureLDA(
      X = train_data$X,
      ICD = train_data$ICD,
      NLP = train_data$NLP,
      HU = train_data$HU,
      filter = train_data$filter,
      prior = "MAP",  # Use internal MAP
      weight = params$weight,
      nEmpty = params$nEmpty,
      alpha = params$alpha,
      beta = params$beta,
      burnin = params$burnin,
      ITER = params$ITER,
      verbose = params$verbose
    )

    # Test using trained phi
    test_results <- sureLDA(
      X = test_data$X,
      ICD = test_data$ICD,
      NLP = test_data$NLP,
      HU = test_data$HU,
      filter = test_data$filter,
      prior = "MAP",  # Use internal MAP
      weight = params$weight,
      nEmpty = params$nEmpty,
      phi = train_results$phi,
      verbose = params$verbose
    )

    cat("sureLDA MAP completed successfully\n")
    list(
      train = train_results,
      test = test_results,
      train_idx = train_idx,
      test_idx = test_idx
    )
  }, error = function(e) {
    warning(paste("Error in sureLDA MAP:", e$message))
    return(NULL)
  })
}

## Run sureLDA with external MAP prior--------------------------------------------------------------------------------------------------
run_surelda_map_p <- function(data, map_results, params = NULL) {
  if (is.null(params)) {
    params <- list(
      weight = "uniform",
      nEmpty = 20,
      alpha = 100,
      beta = 100,
      burnin = 50,
      ITER = 150,
      verbose = TRUE
    )
  }

  train_indices <- map_results$train_idx
  test_indices <- map_results$test_idx

  train_data <- list(
    X = as.matrix(data$X[train_indices,]),
    ICD = as.matrix(data$ICD[train_indices,]),
    NLP = as.matrix(data$NLP[train_indices,]),
    HU = as.vector(data$HU[train_indices]),
    filter = as.matrix(data$filter[train_indices,])
  )

  test_data <- list(
    X = as.matrix(data$X[test_indices,]),
    ICD = as.matrix(data$ICD[test_indices,]),
    NLP = as.matrix(data$NLP[test_indices,]),
    HU = as.vector(data$HU[test_indices]),
    filter = as.matrix(data$filter[test_indices,])
  )

  # Get MAP scores and convert to matrices
  train_prior <- as.matrix(map_results$train$scores)
  test_prior <- as.matrix(map_results$test$scores)

  cat("\nPrior dimensions:")
  cat("\nTrain prior:", dim(train_prior))
  cat("\nTest prior:", dim(test_prior))

  tryCatch({
    train_results <- sureLDA(
      X = train_data$X,
      ICD = train_data$ICD,
      NLP = train_data$NLP,
      HU = train_data$HU,
      filter = train_data$filter,
      prior = train_prior,
      weight = params$weight,
      nEmpty = params$nEmpty,
      alpha = params$alpha,
      beta = params$beta,
      burnin = params$burnin,
      ITER = params$ITER,
      verbose = params$verbose
    )

    test_results <- sureLDA(
      X = test_data$X,
      ICD = test_data$ICD,
      NLP = test_data$NLP,
      HU = test_data$HU,
      filter = test_data$filter,
      prior = test_prior,
      weight = params$weight,
      nEmpty = params$nEmpty,
      phi = train_results$phi,
      verbose = params$verbose
    )

    list(
      train = train_results,
      test = test_results,
      train_idx = train_indices,
      test_idx = test_indices
    )
  }, error = function(e) {
    warning(paste("Error in sureLDA MAP with prior:", e$message))
    print("Error details:")
    print(str(e))
    return(NULL)
  })
}

# Helper function for debugging
debug_matrix <- function(matrix_name, matrix) {
  cat(sprintf("Debug info for %s:\n", matrix_name))
  cat("Dimensions:", dim(matrix), "\n")
  cat("Class:", class(matrix), "\n")
  cat("Is sparse:", inherits(matrix, "sparseMatrix"), "\n")
  cat("Has NAs:", any(is.na(matrix)), "\n")
  if(inherits(matrix, "sparseMatrix")) {
    cat("Non-zero elements:", length(matrix@x), "\n")
  }
}

record_icd_counts <- function(data, test_size = 200, seed = 4747) {
  cat("Recording ICD count values...\n")
  # Split data using outcome vector length
  # split_result <- split_data(data.frame(y = data$outcome), test_size = test_size, seed = seed)

  # Get data for train/test sets
  train_icd <- data$total_ICD_mention[split_result$train_idx]
  test_icd <- data$total_ICD_mention[split_result$test_idx]
  train_outcome <- data$outcome[split_result$train_idx]
  test_outcome <- data$outcome[split_result$test_idx]

  # Instead of fitting a model, just return the values
  cat("ICD count values recorded\n")
  list(
    train = list(icd_counts = train_icd, outcomes = train_outcome),
    test = list(icd_counts = test_icd, outcomes = test_outcome),
    train_idx = split_result$train_idx,
    test_idx = split_result$test_idx
  )
}



## Run ICD logit - Updated implementation
run_icd_logit <- function(dataset, silver_label = "outcome", features = "total_ICD_mention") {
  # Extract train/test split from dataset
  train_idx <- which(dataset$test_set == 0)
  test_idx <- which(dataset$test_set == 1)
  
  train_data <- dataset[train_idx, ]
  test_data <- dataset[test_idx, ]
  
  # Remove rows with missing outcomes
  train_data <- train_data[!is.na(train_data[[silver_label]]), ]
  
  # Fit logistic regression
  formula_str <- paste(silver_label, "~", paste(features, collapse = " + "))
  model <- glm(as.formula(formula_str), data = train_data, family = binomial)
  
  # Predict on both sets
  train_pred <- predict(model, newdata = train_data, type = "response")
  test_pred <- predict(model, newdata = test_data, type = "response")
  
  return(list(
    train = list(probs = train_pred),
    test = list(probs = test_pred),
    train_idx = train_idx,
    test_idx = test_idx,
    model = model
  ))
}

# Note: run_surelda_phenorm and run_surelda_map_p don't need changes 
# since they get train/test indices from their prior results parameters

# Updated calculate_metrics function to properly save optimal_threshold
calculate_metrics <- function(predictions, true_values, true_outcomes,
                              n_total, n_included) {
  n_missing <- n_total - n_included

  predictions <- as.numeric(predictions)
  if(any(is.na(predictions))) {
    warning("NA values found in predictions")
    predictions[is.na(predictions)] <- mean(predictions, na.rm = TRUE)
  }

  true_outcomes <- as.numeric(true_outcomes)

  # Initialize variables
  auc <- precision <- recall <- f1 <- accuracy <- NA
  prob_mse <- prob_mae <- prob_correlation <- optimal_threshold <- NA

  tryCatch({
    roc_obj <- roc(true_outcomes, predictions)
    auc <- as.numeric(auc(roc_obj))
    
    # Get optimal threshold using Youden's J statistic
    coords_result <- coords(roc_obj, "best", best.method = "youden")
    optimal_threshold <- coords_result$threshold
    
    # Create binary predictions using optimal threshold
    binary_preds <- predictions >= optimal_threshold

    # Calculate confusion matrix metrics
    cm <- table(Predicted = binary_preds, Actual = true_outcomes)
    
    # Handle edge cases where confusion matrix might not be 2x2
    if (nrow(cm) == 2 && ncol(cm) == 2) {
      precision <- cm[2,2] / sum(cm[2,])
      recall <- cm[2,2] / sum(cm[,2])
      f1 <- 2 * (precision * recall) / (precision + recall)
      accuracy <- sum(diag(cm)) / sum(cm)
    } else {
      # Handle degenerate cases
      precision <- recall <- f1 <- accuracy <- NA
    }

    # Calculate probability metrics
    prob_mse <- mean((predictions - true_values)^2)
    prob_mae <- mean(abs(predictions - true_values))
    prob_correlation <- cor(predictions, true_values)

  }, error = function(e) {
    warning("Error in metrics calculation: ", e$message)
  })

  data.frame(
    n_total = n_total,
    n_included = n_included,
    n_missing = n_missing,
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

# Fix record_icd_counts function to extract train/test split from data
record_icd_counts <- function(data, test_size = 200, seed = 4747) {
  cat("Recording ICD count values...\n")
  
  # Extract train/test split from data
  train_idx <- which(data$test_set == 0)
  test_idx <- which(data$test_set == 1)

  # Get data for train/test sets
  train_icd <- data$total_ICD_mention[train_idx]
  test_icd <- data$total_ICD_mention[test_idx]
  train_outcome <- data$outcome[train_idx]
  test_outcome <- data$outcome[test_idx]

  # Instead of fitting a model, just return the values
  cat("ICD count values recorded\n")
  list(
    train = list(icd_counts = train_icd, outcomes = train_outcome),
    test = list(icd_counts = test_icd, outcomes = test_outcome),
    train_idx = train_idx,
    test_idx = test_idx
  )
}
