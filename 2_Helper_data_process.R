# 2_Helper_data_process.R----------------------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(Matrix)

##Process Phenorm-------------------------------------------------------------------------------------------------------------------------------------------
process_phenorm <- function(input_data, data_type = c("ori", "ideal", "LDA"), scenario,
                            study_id = "Studyid", valid_label = "Train_Eval_Set",
                            gold_label = "PTYPE_POSITIVE", utilization = "Utiliz",
                            weight = "Sampling_Weight", chart_reviewed = TRUE) {

  data_type <- match.arg(data_type)
  scenario <- match.arg(scenario, c("rare_info", "rare_noninfo", "nonrare_info", "nonrare_noninfo"))

  # Add total_ICD_NLP for ideal and LDA
  if (data_type %in% c("ideal", "LDA")) {
    input_data <- input_data %>%
      mutate(total_ICD_NLP = total_ICD_mention + total_NLP_silver)
  }

  # Set validation labels
  if (!any(grepl(valid_label, names(input_data)))) {
    input_data[[valid_label]] <- rep(0, nrow(input_data))
  }
  if (!is.numeric(input_data %>% pull(valid_label))) {
    valid_label_index <- which(grepl(valid_label, names(input_data), ignore.case = TRUE))
    input_data[[valid_label_index]] <- ifelse(input_data[[valid_label_index]] == "Training", 0, 1)
  }

  # Convert categorical to numeric
  input_data <- input_data %>%
    mutate(across(where(is.character), as.factor)) %>%
    mutate(across(where(is.factor), as.numeric))

  # Rename silver labels
  if (data_type == "ori") {
    input_data <- input_data %>%
      rename(
        total_text_mentions_silver = total_text_mentions,
        total_ICD_NLP_silver = total_ICD_NLP,
        total_ICD_mention_silver = total_ICD_mention,
        total_NLP_silver_silver = total_NLP_silver
      )
  } else {
    input_data <- input_data %>%
      rename(
        total_ICD_NLP_silver = total_ICD_NLP,
        total_ICD_mention_silver = total_ICD_mention,
        total_NLP_silver_silver = total_NLP_silver
      )
  }

  # Process data
  data_processed <- input_data %>%
    mutate(
      PTYPE_POSITIVE = as.numeric(outcome),
      Studyid = as.character(i),
      # Train_Eval_Set = 1,
      Utiliz = 1,
      Sampling_Weight = 1
    ) %>%
    select(-outcome)

  # Get feature names
  data_names <- names(data_processed)
  silver_labels <- data_names[grepl("total", data_names, ignore.case = TRUE) &
                                !grepl("_norm", data_names, ignore.case = TRUE)]
  cui_names <- data_names[grepl("NLP", data_names, ignore.case = TRUE) &
                            grepl("count", data_names, ignore.case = TRUE) &
                            !grepl("case|control", data_names, ignore.case = TRUE)]
  nlp_names <- c(silver_labels, "Utiliz", cui_names)
  structured_data_names <- data_names[!(data_names %in% c(study_id, gold_label,
                                                          valid_label, nlp_names, weight,
                                                          "p_case", "p_control", "i", "note",
                                                          data_names[grepl("probability", data_names, ignore.case = TRUE)]))]

  # Process features
  processed_data <- process_data(
    dataset = data_processed,
    structured_data_names = structured_data_names,
    nlp_data_names = nlp_names,
    study_id = study_id,
    validation_name = valid_label,
    gold_label = gold_label,
    utilization_variable = utilization,
    weight = weight,
    train_on_gold_data = FALSE,
    chart_reviewed = chart_reviewed
  )

  # Filter features
  all_num_unique <- sapply(processed_data$train, function(x) length(unique(x)))
  is_zero_one <- (all_num_unique == 0) | (all_num_unique == 1)
  train <- processed_data$train %>%
    select(!!sym(study_id), all_of(silver_labels), !!sym(weight), !!sym(utilization),
           !!sym(valid_label), (1:ncol(processed_data$train))[!is_zero_one])
  test <- processed_data$test
  outcomes <- processed_data$outcome
  all_data <- processed_data$all

  # Process complete cases
  train_all_cc <- train %>%
    filter(complete.cases(train)) %>%
    select(-!!sym(study_id), -!!sym(valid_label))
  train_screened_cc <- train_all_cc

  test_all_cc <- test %>%
    filter(complete.cases(test)) %>%
    select(-!!sym(study_id), -!!sym(valid_label), -!!sym(gold_label))
  test_screened_cc <- test_all_cc

  # Log transform
  analysis_data <- list(
    "train" = apply_log_transformation(train_screened_cc,
                                       names(train_screened_cc)[!grepl(weight, names(train_screened_cc))],
                                       utilization),
    "test" = apply_log_transformation(test_screened_cc,
                                      names(test_screened_cc)[!grepl(weight, names(test_screened_cc))],
                                      utilization),
    "outcomes" = outcomes,
    "train_all" = apply_log_transformation(train_all_cc,
                                           names(train_all_cc)[!grepl(weight, names(train_all_cc))],
                                           utilization),
    "test_all" = apply_log_transformation(test_all_cc,
                                          names(test_all_cc)[!grepl(weight, names(test_all_cc))],
                                          utilization),
    "all" = apply_log_transformation(all_data,
                                     names(all_data)[!(names(all_data) %in%
                                                         c(gold_label, valid_label, study_id, weight))],
                                     utilization),
    "utilization_variable" = utilization,
    "silver_labels" = silver_labels
  )

  return(analysis_data)
}

## Updated process_map - now uses pre-generated note variable
process_map <- function(df, data_type = c("ori", "ideal", "LDA"), scenario) {
  data_type <- match.arg(data_type)

  # Process counts
  processed <- df %>%
    mutate(
      total_ICD_mention = round(total_ICD_mention),
      total_NLP_silver = round(total_NLP_silver)
      # Remove note generation - now using pre-generated note variable
    )

  # Create matrix format
  mat <- Matrix(as.matrix(processed[, c("total_ICD_mention", "total_NLP_silver")]), sparse=TRUE)
  colnames(mat) <- c("ICD", "NLP")

  # Use the pre-generated note variable
  note <- Matrix(df$note, sparse=TRUE)
  colnames(note) <- "note"

  return(list(mat=mat, note=note, test_set=df$test_set))
}

## Process SureLDA - Updated to include test_set
process_surelda <- function(df, data_type = c("ori", "ideal", "LDA"), scenario) {
  data_type <- match.arg(data_type)

  # Create basic matrices
  processed <- list()

  # Generate feature matrix based on data type
  if (data_type == "ideal") {
    processed$X <- Matrix(cbind(
      df$total_ICD_mention,
      df$total_NLP_silver
    ), sparse=TRUE)
    colnames(processed$X) <- c("ICD", "NLP")
  } else if (data_type == "LDA") {
    nlp_cols <- grep("^NLP_[0-9]+$", names(df), value=TRUE)
    processed$X <- Matrix(as.matrix(cbind(
      df$total_ICD_mention,
      df$total_NLP_silver,
      df[, nlp_cols]
    )), sparse=TRUE)
    colnames(processed$X) <- c("ICD", "NLP", nlp_cols)
  } else {
    nlp_count_cols <- grep("^NLP_[0-9]+_count$", names(df), value=TRUE)
    selected_count_cols <- grep("^Selected_NLP_[0-9]+_count$", names(df), value=TRUE)
    processed$X <- Matrix(as.matrix(cbind(
      df$total_ICD_mention,
      df$total_NLP_silver,
      df[, c(nlp_count_cols, selected_count_cols)]
    )), sparse=TRUE)
    colnames(processed$X) <- c("ICD", "NLP", nlp_count_cols, selected_count_cols)
  }

  # Create individual matrices
  processed$ICD <- Matrix(df$total_ICD_mention, sparse=TRUE)
  processed$NLP <- Matrix(df$total_NLP_silver, sparse=TRUE)
  processed$filter <- Matrix((df$total_ICD_mention > 0) * 1, sparse=TRUE)
  processed$HU <- Matrix(
    if (data_type == "ori") {
      df$n_i_num
    } else if (data_type == "ideal") {
      df$note
    } else {
      df$H
    },
    sparse = TRUE
  )
  if("outcome" %in% names(df)) {
    processed$outcome <- df$outcome
  }
  
  # Add test_set information
  processed$test_set <- df$test_set

  return(processed)
}

## Updated process_surelda_map - now uses pre-generated note variable
process_surelda_map <- function(df, data_type = c("ori", "ideal", "LDA"), scenario) {
  data_type <- match.arg(data_type)

  # Round counts
  df <- df %>%
    mutate(
      total_ICD_mention = round(total_ICD_mention),
      total_NLP_silver = round(total_NLP_silver)
    )

  processed <- list()

  # Create feature matrix based on data type
  if (data_type == "ideal") {
    processed$X <- Matrix(cbind(
      df$total_ICD_mention,
      df$total_NLP_silver
    ), sparse=TRUE)
    colnames(processed$X) <- c("ICD", "NLP")
  } else {
    features <- if(data_type == "LDA") {
      nlp_cols <- grep("^NLP_[0-9]+$", names(df), value=TRUE)
      round(as.matrix(df[, nlp_cols]))
    } else {
      nlp_count_cols <- grep("^NLP_[0-9]+_count$", names(df), value=TRUE)
      selected_count_cols <- grep("^Selected_NLP_[0-9]+_count$", names(df), value=TRUE)
      round(as.matrix(df[, c(nlp_count_cols, selected_count_cols)]))
    }

    processed$X <- Matrix(cbind(
      df$total_ICD_mention,
      df$total_NLP_silver,
      features
    ), sparse=TRUE)
  }

  # Create individual matrices
  processed$ICD <- Matrix(df$total_ICD_mention, sparse=TRUE)
  processed$NLP <- Matrix(df$total_NLP_silver, sparse=TRUE)
  processed$filter <- Matrix((df$total_ICD_mention > 0) * 1, sparse=TRUE)
  
  # Use the pre-generated note variable instead of generating new one
  processed$HU <- Matrix(df$note, sparse=TRUE)

  if("outcome" %in% names(df)) {
    processed$outcome <- df$outcome
  }
  
  # Add test_set information
  processed$test_set <- df$test_set

  return(processed)
}

## Helper functions------------------------------------------------------------------------------------------------------------------------------------------
process_structured_data <- function(data, vars_to_process, values) {
  # is_chr <-
  is_chr_bin <- apply(data, 2, function(x) length(unique(x)) == 2 & !is.numeric(x))
  bin_names <- names(data)[is_chr_bin]
  for (i in seq_len(length(bin_names))) {
    var <- bin_names[i]
    val <- values[i]
    indx <- which(grepl(var, names(data), ignore.case = TRUE))
    data[[indx]] <- ifelse(data[[indx]] == val, 1, 0)
  }
  return(data)
}

process_data <- function(dataset = NULL, structured_data_names = "AGE",
                         nlp_data_names = "C", study_id = "STUDYID",
                         validation_name = "GOLD_STANDARD_VALIDATION",
                         gold_label = "AP_GOLD_LABEL",
                         utilization_variable = "Utiliz",
                         weight = "weight",
                         train_on_gold_data = FALSE,
                         chart_reviewed = TRUE) {
  if (!any(grepl(utilization_variable, names(dataset)))) {
    dataset[[utilization_variable]] <- 1
  }
  if (!any(grepl(weight, names(dataset)))) {
    dataset[[weight]] <- 1
  }
  if (chart_reviewed) {
    all_data <- dplyr::select(
      dataset, !!c(matches(study_id), matches(validation_name),
                   matches(gold_label), matches(weight),
                   matches(paste0("^", structured_data_names, "$")),
                   matches(unique(c(nlp_data_names, utilization_variable))))
    )
    outcome_indx <- which(grepl(gold_label, names(all_data), ignore.case = TRUE))
    valid_indx <- which(grepl(validation_name, names(all_data), ignore.case = TRUE))
    outcomes <- all_data[[outcome_indx]]
    if (!is.numeric(outcomes)) {
      if (any(grepl("yes", outcomes, ignore.case = TRUE))) {
        outcomes <- ifelse(grepl("yes", outcomes, ignore.case = TRUE), 1, 0)
        outcomes[all_data[[valid_indx]] == 0] <- NA
      } else {
        stop("Outcome is not in a recognized format (numeric, binary, or 'yes'/'no'). Please use an outcome variable that is in one of these formats.")
      }
      all_data[[outcome_indx]] <- outcomes
    }
    if (train_on_gold_data) {
      train <- dplyr::select(all_data, -!!matches(gold_label))
    } else {
      train <- dplyr::select(dplyr::filter(all_data, !!rlang::sym(names(all_data)[valid_indx]) == 0),
                             -!!matches(gold_label))
    }
    test <- dplyr::filter(all_data, !!rlang::sym(names(all_data)[valid_indx]) == 1)
  } else {
    all_data <- dplyr::select(
      dataset, !!c(matches(study_id), matches(validation_name), matches(weight),
                   matches(paste0("^", structured_data_names, "$")),
                   matches(unique(c(nlp_data_names, utilization_variable))))
    )
    train <- all_data
    test <- all_data
  }
  train_cc <- train[complete.cases(train), ]
  test_cc <- test[complete.cases(test), ]
  if (chart_reviewed) {
    outcome_indx <- which(grepl(gold_label, names(test_cc), ignore.case = TRUE))
    outcomes <- test_cc[[outcome_indx]]
  } else {
    outcomes <- rep(NA, nrow(test_cc))
  }
  return(list("outcome" = outcomes, "train" = train_cc, "test" = test_cc,
              "all" = all_data))
}

apply_log_transformation <- function(dataset = NULL, varnames = NULL,
                                     utilization_var = NULL) {
  log_transformed_data <- dataset
  log_transformed_data[, varnames] <- log(dataset[, varnames] + 1)
  if (all(dataset[[utilization_var]] == 1)) {
    log_transformed_data[[utilization_var]] <- 1
  }
  return(log_transformed_data)
}

process_surelda_multi_silver <- function(df, data_type = c("ori", "ideal", "LDA"), scenario, n_silver = 4) {
  data_type <- match.arg(data_type)
  
  processed <- process_surelda(df, data_type, scenario)
  
  data_names <- names(df)
  silver_label_vars <- data_names[grepl("total", data_names, ignore.case = TRUE) &
                                    !grepl("_norm", data_names, ignore.case = TRUE)]
  
  if (length(silver_label_vars) < 2) {
    stop("Need at least 2 silver label variables, found: ", length(silver_label_vars))
  }
  
  if (n_silver > length(silver_label_vars)) {
    warning("Requested n_silver (", n_silver, ") > available silver labels (", 
            length(silver_label_vars), "). Using all available.")
    n_silver <- length(silver_label_vars)
  }
  
  selected_vars <- silver_label_vars[1:n_silver]
  
  silver_labels <- list()
  n_patients <- nrow(df)
  
  for (i in 1:length(selected_vars)) {
    var_name <- selected_vars[i]
    clean_name <- gsub("^total_", "", var_name, ignore.case = TRUE)
    clean_name <- gsub("_silver$", "", clean_name, ignore.case = TRUE)
    
    silver_labels[[clean_name]] <- as.numeric(df[[var_name]])
    
    if (length(silver_labels[[clean_name]]) != n_patients) {
      stop(paste("Silver label", clean_name, "has length", length(silver_labels[[clean_name]]), 
                 "but expected", n_patients))
    }
  }
  
  lengths <- sapply(silver_labels, length)
  if (!all(lengths == n_patients)) {
    stop("All silver labels must have the same length as the number of patients")
  }
  
  processed$silver_labels <- silver_labels
  processed$n_silver <- length(silver_labels)
  processed$silver_label_vars <- selected_vars
  
  return(processed)
}