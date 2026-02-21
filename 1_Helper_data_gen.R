### Data generation
## Updated create_simulation_data (for "ideal" data type)
create_simulation_data <- function(
    scenario = c("rare_info", "rare_noninfo", "nonrare_info", "nonrare_noninfo"),
    n = 10000, test_size = 200,
    rare_prob = list(mean = 0.05, sd = 0.1),
    nonrare_prob = list(mean = 0.4, sd = 0.1),
    bimodal_counts = list(
      case = list(lambda = 8),
      control = list(lambda = 2)
    ),
    unimodal_counts = list(
      case = list(lambda = 5),
      control = list(lambda = 5)
    ),
    count2_multiplier = 1.5) {
  require(dplyr)
  scenario <- match.arg(scenario)
  params <- switch(scenario,
                   rare_info = list(prob = rare_prob, counts = bimodal_counts),
                   rare_noninfo = list(prob = rare_prob, counts = unimodal_counts),
                   nonrare_info = list(prob = nonrare_prob, counts = bimodal_counts),
                   nonrare_noninfo = list(prob = nonrare_prob, counts = unimodal_counts))
  initial_probability <- pmax(pmin(rnorm(n, params$prob$mean, params$prob$sd), 1), 0)
  initial_data <- data.frame(
    i = 1:n,
    initial_probability = initial_probability
  ) %>%
    rowwise() %>%
    mutate(
      outcome = rbinom(1, 1, initial_probability),
      note = 1 + rnorm(1, mean=0, sd=0.001)  # Generate noisy note variable here
    ) %>%
    mutate(
      total_ICD_mention = if(outcome == 1) {
        rpois(1, params$counts$case$lambda)
      } else {
        rpois(1, params$counts$control$lambda)
      },
      total_NLP_silver = if(outcome == 1) {
        rpois(1, params$counts$case$lambda * count2_multiplier)
      } else {
        rpois(1, params$counts$control$lambda * count2_multiplier)
      }
    ) %>%
    mutate(
      p_case = initial_probability *
        dpois(total_ICD_mention, params$counts$case$lambda) *
        dpois(total_NLP_silver, params$counts$case$lambda * count2_multiplier),
      p_control = (1 - initial_probability) *
        dpois(total_ICD_mention, params$counts$control$lambda) *
        dpois(total_NLP_silver, params$counts$control$lambda * count2_multiplier),
      probability = if_else(p_case + p_control == 0,
                            initial_probability,
                            p_case / (p_case + p_control))
    ) %>%
    ungroup()
  # get train/test split
  train_test <- split_data(initial_data, test_size = test_size)
  test_set <- rep(0, nrow(initial_data))
  test_set[train_test$test_idx] <- 1
  return(initial_data %>% mutate(test_set = test_set))
}

## Updated create_LDA_data (for "LDA" data type)
create_LDA_data <- function(n = 10000, test_size = 200,
                            scenario = c("rare_info", "rare_noninfo",
                                         "nonrare_info", "nonrare_noninfo"),
                            prob_rare = list(mean = -0.83, sd = sqrt(0.25)),
                            prob_nonrare = list(mean = -0.095, sd = sqrt(0.25)),
                            beta_info = list(
                              beta_plus = list(ICD = 5.5, NLP = 2.1),
                              beta_minus = list(ICD = 0.71, NLP = 0.81)
                            ),
                            beta_noninfo = list(
                              beta_plus = list(ICD = 3.0, NLP = 1.0),
                              beta_minus = list(ICD = 3.0, NLP = 1.0)
                            ),
                            output_dir = NULL) {

  scenario <- match.arg(scenario)

  params <- switch(scenario,
                   rare_info = list(prob = prob_rare, beta = beta_info),
                   rare_noninfo = list(prob = prob_rare, beta = beta_noninfo),
                   nonrare_info = list(prob = prob_nonrare, beta = beta_info),
                   nonrare_noninfo = list(prob = prob_nonrare, beta = beta_noninfo))

  # Generate outcome
  Y_continuous <- rnorm(n, mean = params$prob$mean, sd = params$prob$sd)
  Y_binary <- ifelse(Y_continuous > 0, 1, -1)

  # Generate visits
  H <- pmax(rpois(n, lambda = 2), 1)
  
  # Generate noisy note variable for MAP consistency
  note_var <- 1 + rnorm(n, mean=0, sd=0.001)
  
  demographics <- data.frame(
    i = 1:n, 
    H = H,
    note = note_var  # Add note variable here
  )

  # Generate ICD mentions
  ICD_count <- sapply(1:n, function(i) {
    beta_ICD <- if(Y_binary[i] == 1) params$beta$beta_plus$ICD else params$beta$beta_minus$ICD
    count <- (rgamma(1, beta_ICD, 1) * H[i]^0.3 - 1)
    if(Y_binary[i] == 1) max(count, 1) else max(count, 0)
  })

  # Generate NLP mentions
  NLP_count <- sapply(1:n, function(i) {
    beta_NLP <- if(Y_binary[i] == 1) params$beta$beta_plus$NLP else params$beta$beta_minus$NLP
    count <- (rgamma(1, beta_NLP, 1) * H[i]^0.25 - 1)
    if(Y_binary[i] == 1) max(count, 1) else max(count, 0)
  })

  # Create result dataframe
  result <- demographics %>%
    mutate(
      total_ICD_mention = ICD_count,
      total_NLP_silver = NLP_count
    )

  # Generate additional NLP features
  for(j in 1:150) {
    result[[paste0("NLP_", j)]] <- sapply(1:n, function(i) {
      b_j <- if(Y_binary[i] == 1) 1.25 else 0.9
      count <- (rgamma(1, b_j, 1) * H[i]^0.2 - 1)
      max(count, 0)
    })
  }

  # Calculate probabilities
  result <- result %>%
    mutate(
      outcome = ifelse(Y_binary == -1, 0, 1),
      probability = sapply(1:n, function(i) {
        prob <- 1/(1 + exp(-Y_continuous[i]))

        log_case <- dgamma(ICD_count[i], params$beta$beta_plus$ICD, 1, log=TRUE) +
          dgamma(NLP_count[i], params$beta$beta_plus$NLP, 1, log=TRUE)

        log_control <- dgamma(ICD_count[i], params$beta$beta_minus$ICD, 1, log=TRUE) +
          dgamma(NLP_count[i], params$beta$beta_minus$NLP, 1, log=TRUE)

        log_num <- log(prob) + log_case
        log_denom <- log(prob * exp(log_case) + (1-prob) * exp(log_control))

        exp(log_num - log_denom)
      })
    )

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    filename <- paste0("input_data_LDA_", scenario, ".csv")
    write.csv(result, file.path(output_dir, filename), row.names = FALSE)
  }
  # get train/test split
  train_test <- split_data(result, test_size = test_size)
  test_set <- rep(0, nrow(result))
  test_set[train_test$test_idx] <- 1
  return(result %>% mutate(test_set = test_set))
}

## Updated create_complex_data (for "ori" data type)
create_complex_data <- function(n = 10000,
                                test_size = 200,
                                scenario = c("rare_info", "rare_noninfo",
                                             "nonrare_info", "nonrare_noninfo"),
                                rare_disease = list(
                                  beta_Y = c(-5.7, 0.05, 0.03, 0.1, rep(0.2, 7)),
                                  lambda = 5
                                ),
                                nonrare_disease = list(
                                  beta_Y = c(-2.8, 0.05, 0.03, 0.1, rep(0.2, 7)),
                                  lambda = 5
                                ),
                                info_labels = list(
                                  offset_case = 0.4,
                                  offset_control = 0.1,
                                  nlp_case_offset = 0.3,
                                  nlp_control_offset = 0.1,
                                  icd_case_offset = 0.3,
                                  icd_control_offset = 0.1,
                                  Coef_T = c(-2, 0.03, -0.2, 0.1, 0.5, -0.5, 0.8, -0.5, 0.00, 0.2, 0.5),
                                  Coef_S = c(-2.5, 0.02, -0.1, 0.5, -0.02, 0.00, -0.3, 0.03, 0.00, 0.7, -0.3)
                                ),
                                noninfo_labels = list(
                                  offset_case = 0.1,
                                  offset_control = 0.1,
                                  nlp_case_offset = 0.1,
                                  nlp_control_offset = 0.1,
                                  icd_case_offset = 0.1,
                                  icd_control_offset = 0.1,
                                  Coef_T = c(-2, 0.03, -0.2, 0.1, 0.5, -0.5, 0.8, -0.5, 0.00, 0.2, 0.5),
                                  Coef_S = c(-2.5, 0.02, -0.1, 0.5, -0.02, 0.00, -0.3, 0.03, 0.00, 0.7, -0.3)
                                ),
                                base_zero_prob = 0.90,
                                variation_range = 0.1,
                                lambda_rate = 1/10,
                                output_dir = NULL) {

  scenario <- match.arg(scenario)

  params <- switch(scenario,
                   rare_info = list(disease_params = rare_disease, silver_params = info_labels),
                   rare_noninfo = list(disease_params = rare_disease, silver_params = noninfo_labels),
                   nonrare_info = list(disease_params = nonrare_disease, silver_params = info_labels),
                   nonrare_noninfo = list(disease_params = nonrare_disease, silver_params = noninfo_labels))

  adjust_prob <- function(prob) {
    pmax(pmin(replace(prob, is.na(prob), 0.001), 1), 0.001)
  }

  # Demographics
  race_levels <- c("Asian", "Black or African American", "Native American or Alaska Native",
                   "Native Hawaiian or Pacific Islander", "Other", "Unknown", "White")
  race_probs <- c(0.13, 0.073, 0.017, 0.015, 0.034, 0.12, 0.61)

  initial_data <- data.frame(
    i = 1:n,
    age = rnorm(n, mean=40.7, sd=22.6),
    sex = rbinom(n, 1, prob=0.51),
    race = sample(race_levels, size=n, replace=TRUE, prob=race_probs),
    exposure = rbinom(n, 1, prob=0.30),
    n_i_num = rpois(n, lambda=params$disease_params$lambda) + 1,
    note = 1 + rnorm(n, mean=0, sd=0.001)  # Generate noisy note variable here
  ) %>%
    filter(age > 0 & n_i_num > 0)

  # Generate outcome
  X_race <- model.matrix(~ race - 1, data=initial_data)
  X_demo <- model.matrix(~ age + sex + exposure, data=initial_data)
  X <- cbind(X_demo, X_race)

  logit_prob_Y <- as.vector(X %*% params$disease_params$beta_Y)
  initial_data$P_initial <- 1 / (1 + exp(-logit_prob_Y))
  initial_data$Outcome <- rbinom(nrow(initial_data), 1, prob=initial_data$P_initial)

  # Text probabilities
  logit_P_T <- as.vector(X %*% params$silver_params$Coef_T)
  initial_data$P_Ti <- 1 / (1 + exp(-logit_P_T))

  logit_P_S <- as.vector(X %*% params$silver_params$Coef_S)
  initial_data$P_Si <- 1 / (1 + exp(-logit_P_S))

  # Add offsets
  initial_data <- initial_data %>%
    rowwise() %>%
    mutate(
      P_T0 = adjust_prob(P_Ti + params$silver_params$offset_control),
      P_T1 = adjust_prob(P_Ti + params$silver_params$offset_case),
      P_S0 = adjust_prob(P_Si + params$silver_params$offset_control),
      P_S1 = adjust_prob(P_Si + params$silver_params$offset_case),
      P_T = adjust_prob(ifelse(Outcome == 1,
                               P_Ti + params$silver_params$offset_case,
                               P_Ti + params$silver_params$offset_control)),
      P_S = adjust_prob(ifelse(Outcome == 1,
                               P_Si + params$silver_params$offset_case,
                               P_Si + params$silver_params$offset_control))
    ) %>%
    ungroup()

  initial_data$max_counts <- sample(1:15, nrow(initial_data), replace=TRUE)

  # Expand data
  expanded_data <- initial_data %>%
    slice(rep(1:n(), initial_data$n_i_num))
  expanded_data$n_i <- unlist(lapply(initial_data$n_i_num, seq_len))

  # Generate text strings
  expanded_data <- expanded_data %>%
    mutate(
      Text_string_T = rbinom(n(), max_counts, prob=P_T),
      Text_string_S = rbinom(n(), max_counts, prob=P_S),
      Text_string = Text_string_T + Text_string_S
    )

  # NLP Variables
  lambda_values <- rexp(150, rate=lambda_rate)
  zero_probs <- runif(150, min=base_zero_prob - variation_range,
                      max=base_zero_prob + variation_range)

  expanded_data <- cbind(expanded_data,
                         matrix(NA, nrow=nrow(expanded_data), ncol=150))
  colnames(expanded_data)[(ncol(expanded_data)-149):ncol(expanded_data)] <- paste0("NLP_", 1:150)

  for(j in 1:150) {
    expanded_data[[paste0("NLP_", j)]] <- ifelse(runif(nrow(expanded_data)) < zero_probs[j],
                                                 0,
                                                 ceiling(rexp(nrow(expanded_data), rate=lambda_values[j])))
  }

  # NLP silver mentions
  expanded_data <- expanded_data %>%
    mutate(
      P_NLP_silver = adjust_prob((P_T + P_S)/2),
      P_NLP_silver_1 = adjust_prob((P_T1 + P_S1)/2),
      P_NLP_silver_0 = adjust_prob((P_T0 + P_S0)/2),
      NLP_silver = rbinom(n(), 1, prob=P_NLP_silver),
      NLP_silver_case = rbinom(n(), 1, prob=P_NLP_silver_1),
      NLP_silver_control = rbinom(n(), 1, prob=P_NLP_silver_0),
      NLP_silver_count = ifelse(NLP_silver == 1,
                                ceiling(rexp(n(), rate=1/P_NLP_silver)),
                                0)
    )

  # Selected NLP features
  generate_coefficients <- function() {
    c(runif(1, -4,-2), runif(1, -0.05, 0.05), runif(1, -0.1, 0.1),
      runif(1, -0.5, 0.5), runif(1, -0.5, 0.5), runif(6, -0.3, 0.3),
      runif(1, 0, 0.2))
  }

  coefficients_list <- replicate(10, generate_coefficients(), simplify=FALSE)
  coef_nlp <- do.call(rbind, coefficients_list)

  for(k in 1:10) {
    selected_nlp_col <- paste0("Selected_NLP_", k)
    selected_nlp_case_col <- paste0("Selected_NLP_case", k)
    selected_nlp_control_col <- paste0("Selected_NLP_control", k)

    demographic_component <- as.vector(X %*% coef_nlp[k, 2:12])
    Pj_base <- 1 / (1 + exp(-demographic_component))

    Pj <- adjust_prob(ifelse(expanded_data$Outcome == 1,
                             Pj_base + params$silver_params$nlp_case_offset,
                             Pj_base - params$silver_params$nlp_control_offset))
    pj_0 <- adjust_prob(Pj_base - params$silver_params$nlp_control_offset)
    pj_1 <- adjust_prob(Pj_base + params$silver_params$nlp_case_offset)

    expanded_data[[selected_nlp_col]] <- rbinom(nrow(expanded_data), 1, prob=Pj)
    expanded_data[[selected_nlp_case_col]] <- rbinom(nrow(expanded_data), 1, prob=pj_1)
    expanded_data[[selected_nlp_control_col]] <- rbinom(nrow(expanded_data), 1, prob=pj_0)

    expanded_data[[paste0(selected_nlp_col, "_count")]] <-
      ifelse(expanded_data[[selected_nlp_col]] == 1,
             ceiling(rpois(nrow(expanded_data), lambda=Pj)),
             0)
    expanded_data[[paste0(selected_nlp_case_col, "_count")]] <-
      ifelse(expanded_data[[selected_nlp_case_col]] == 1,
             ceiling(rpois(nrow(expanded_data), lambda=pj_1)),
             0)
    expanded_data[[paste0(selected_nlp_control_col, "_count")]] <-
      ifelse(expanded_data[[selected_nlp_control_col]] == 1,
             ceiling(rpois(nrow(expanded_data), lambda=pj_0)),
             0)
  }

  # ICD mentions
  Coef_ICD <- c(-0.5, runif(10,-0.2,0.2), 0.3)

  X_ICD <- cbind(1, expanded_data %>%
                   select(starts_with("Selected_NLP_") &
                            matches("count") &
                            !matches("case|control")),
                 expanded_data$NLP_silver)
  X_ICD_case <- cbind(1, expanded_data %>%
                        select(starts_with("Selected_NLP_") &
                                 matches("count") &
                                 matches("case")),
                      expanded_data$NLP_silver_case)
  X_ICD_control <- cbind(1, expanded_data %>%
                           select(starts_with("Selected_NLP_") &
                                    matches("count") &
                                    matches("control")),
                         expanded_data$NLP_silver_control)

  logit_ICD <- as.matrix(X_ICD) %*% Coef_ICD
  logit_ICD_case <- as.matrix(X_ICD_case) %*% Coef_ICD
  logit_ICD_control <- as.matrix(X_ICD_control) %*% Coef_ICD

  P_ICD <- (1 / (1 + exp(-logit_ICD)))
  P_ICD <- adjust_prob(ifelse(expanded_data$Outcome == 1,
                              P_ICD + params$silver_params$icd_case_offset,
                              P_ICD - params$silver_params$icd_control_offset))
  P_ICD_case <- adjust_prob((1 / (1 + exp(-logit_ICD_case))) +
                              params$silver_params$icd_case_offset)
  P_ICD_control <- adjust_prob((1 / (1 + exp(-logit_ICD_control))) -
                                 params$silver_params$icd_control_offset)

  expanded_data <- expanded_data %>%
    mutate(
      P_ICD = P_ICD,
      P_ICD_case = P_ICD_case,
      P_ICD_control = P_ICD_control,
      ICD_mention = rbinom(n(), 1, prob=P_ICD),
      ICD_NLP = ICD_mention + NLP_silver
    )

  # Summarize - make sure to include note in the summarization
  summarized_data <- expanded_data %>%
    group_by(i) %>%
    summarise(
      n_i_num = first(n_i_num),
      note = first(note),  # Keep the note variable
      total_NLP_silver = sum(NLP_silver),
      total_text_mentions = sum(Text_string),
      total_ICD_NLP = sum(ICD_NLP),
      total_ICD_mention = sum(ICD_mention),
      age = first(age),
      sex = first(sex),
      race = first(race),
      exposure = first(exposure),
      outcome = first(Outcome),
      P_initial = first(P_initial),
      across(starts_with("NLP_") & !contains("silver"), ~sum(., na.rm=TRUE)),
      across(starts_with("selected_"), ~sum(., na.rm=TRUE)),
      NLP_silver_count = sum(NLP_silver_count),
      Text_string_T = sum(Text_string_T),
      Text_string_S = sum(Text_string_S),
      max_counts = sum(max_counts),
      across(starts_with("P_"), ~mean(., na.rm=TRUE)),
      .groups='drop'
    ) %>%
    rename_with(~paste0("mean_", .), starts_with("P_"))

  # Helper functions for probability calculations
  calculate_PZ_XY <- function(total_text_mentions, n_i_num, mean_P_T, mean_P_S,
                              total_NLP_silver, mean_P_NLP_silver, total_ICD_mention,
                              mean_P_ICD, total_ICD_NLP, max_counts) {
    P_total_text_string <- (mean_P_T + mean_P_S)/2 + 0.001
    P_total_NLP_silver <- mean_P_NLP_silver + 0.001
    P_total_ICD_mention <- mean_P_ICD + 0.001
    P_total_ICD_NLP <- (mean_P_ICD + mean_P_NLP_silver)/2 + 0.001
    PZ_XY <- P_total_text_string * P_total_NLP_silver * P_total_ICD_NLP * P_total_ICD_mention + 0.001
    return(PZ_XY)
  }

  calculate_PX <- function(age, sex, n_i_num, exposure, race) {
    race_levels <- c("Asian", "Black or African American", "Native American or Alaska Native",
                     "Native Hawaiian or Pacific Islander", "Other", "Unknown", "White")
    race_probs <- c(0.13, 0.073, 0.017, 0.015, 0.034, 0.12, 0.61)
    mu_age <- 40.7
    sd_age <- 22.6
    P_age <- dnorm(age, mean=mu_age, sd=sd_age) + 0.001
    P_sex <- 0.51 + 0.001
    P_exposure <- 0.30 + 0.001
    P_race <- race_probs[which(race_levels == race)] + 0.001
    P_X <- P_age * P_sex * P_exposure * P_race + 0.001
    return(P_X)
  }

  # Calculate probabilities
  summarized_data <- summarized_data %>%
    rowwise() %>%
    mutate(
      PZ_XY = calculate_PZ_XY(total_text_mentions, n_i_num, mean_P_T, mean_P_S,
                              total_NLP_silver, mean_P_NLP_silver, total_ICD_mention,
                              mean_P_ICD, total_ICD_NLP, max_counts),
      P_X = calculate_PX(age, sex, n_i_num, exposure, race),
      PZ_XY_0 = calculate_PZ_XY(total_text_mentions, n_i_num,
                                mean_P_T0, mean_P_S0,
                                total_NLP_silver, mean_P_NLP_silver_0,
                                total_ICD_mention, mean_P_ICD_control,
                                total_ICD_NLP, max_counts),
      PZ_XY_1 = calculate_PZ_XY(total_text_mentions, n_i_num,
                                mean_P_T1, mean_P_S1,
                                total_NLP_silver, mean_P_NLP_silver_1,
                                total_ICD_mention, mean_P_ICD_case,
                                total_ICD_NLP, max_counts)
    ) %>%
    ungroup()

  # Generate marginal PY
  n_samples <- 10000
  outcomes <- numeric(n_samples)
  for (i in 1:n_samples) {
    sampled_row <- summarized_data[sample(nrow(summarized_data), 1), ]
    outcomes[i] <- sampled_row$outcome
  }
  P_Y <- mean(outcomes)
  summarized_data$P_Y <- P_Y

  # Calculate final probabilities
  calculate_final_prob <- function(mean_P_initial, PZ_XY, PZ_XY_0, PZ_XY_1, P_Y, P_X) {
    numerator <- (PZ_XY_1 * mean_P_initial * P_X)
    denominator <- (PZ_XY_0 * (1-P_Y) * P_X) + (PZ_XY_1 * P_Y * P_X)
    P_real <- numerator/denominator
    return(P_real)
  }

  summarized_data <- summarized_data %>%
    rowwise() %>%
    mutate(
      P_real = calculate_final_prob(mean_P_initial, PZ_XY, PZ_XY_0, PZ_XY_1, P_Y, P_X)
    ) %>%
    ungroup()

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

    filename <- switch(scenario,
                       "rare_info" = "input_data_ori_rare_info.csv",
                       "rare_noninfo" = "input_data_ori_rare_noninfo.csv",
                       "nonrare_info" = "input_data_ori_nonrare_info.csv",
                       "nonrare_noninfo" = "input_data_ori_nonrare_noninfo.csv")

    write.csv(summarized_data, file=file.path(output_dir, filename), row.names=FALSE)
  }
  #rescale P_real
  summarized_data$probability <- (summarized_data$P_real - min(summarized_data$P_real)) /
    (max(summarized_data$P_real) - min(summarized_data$P_real))
  # get train/test split
  train_test <- split_data(summarized_data, test_size = test_size)
  test_set <- rep(0, nrow(summarized_data))
  test_set[train_test$test_idx] <- 1
  return(summarized_data %>% mutate(test_set = test_set))
}

# Updated processing functions to use the pre-generated note variable

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
