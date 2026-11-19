# ============================================================
# FOUR-STATE MARKOV CHAIN FOR DAILY WEATHER
# Real Data: Wilmington International Airport
# NOAA Daily Weather Data, 2016-2025
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. Read and prepare the data
# ------------------------------------------------------------

weather <- read.csv(file.choose())

weather <- weather %>%
  mutate(
    DATE = as.Date(DATE)
  ) %>%
  arrange(DATE)

head(weather)
dim(weather)

# ------------------------------------------------------------
# 2. Split data into training and testing periods
# Training: 2016-2023
# Testing: 2024-2025
# ------------------------------------------------------------

weather$year <-
  as.integer(format(weather$DATE, "%Y"))

train <- weather %>%
  filter(year <= 2023)

test <- weather %>%
  filter(year >= 2024)

cat(
  "Training observations:",
  nrow(train),
  "\n"
)

cat(
  "Testing observations:",
  nrow(test),
  "\n"
)

# ------------------------------------------------------------
# 3. Select preliminary Heavy Rain threshold
# ------------------------------------------------------------

heavy_threshold <- as.numeric(
  quantile(
    train$PRCP[train$PRCP > 0],
    probs = 0.80,
    na.rm = TRUE
  )
)

cat(
  "Heavy Rain threshold =",
  round(heavy_threshold, 3),
  "inches\n"
)

cat(
  "Heavy Rain threshold =",
  round(heavy_threshold * 25.4, 2),
  "mm\n"
)

# ------------------------------------------------------------
# 4. Define four weather states
# ------------------------------------------------------------

states <- c(
  "Dry",
  "Rain",
  "Heavy Rain",
  "Thunder-Reported Day"
)

classify_weather <-
  function(prcp, thunder, threshold) {

    state <-
      rep(NA_character_, length(prcp))

    storm_day <-
      !is.na(thunder) & thunder == 1

    state[storm_day] <- "Thunder-Reported Day"

    nonstorm <- !storm_day

    state[
      nonstorm &
      !is.na(prcp) &
      prcp == 0
    ] <- "Dry"

    state[
      nonstorm &
      !is.na(prcp) &
      prcp > 0 &
      prcp < threshold
    ] <- "Rain"

    state[
      nonstorm &
      !is.na(prcp) &
      prcp >= threshold
    ] <- "Heavy Rain"

    factor(
      state,
      levels = states
    )
  }

train$State <-
  classify_weather(
    train$PRCP,
    train$WT03,
    heavy_threshold
  )

test$State <-
  classify_weather(
    test$PRCP,
    test$WT03,
    heavy_threshold
  )

# ------------------------------------------------------------
# 5. State distributions
# ------------------------------------------------------------

table(train$State)

round(
  prop.table(
    table(train$State)
  ),
  3
)

table(test$State)

# ------------------------------------------------------------
# 6. Transition count matrix
# ------------------------------------------------------------

current_state <-
  head(train$State, -1)

next_state <-
  tail(train$State, -1)

current_date <-
  head(train$DATE, -1)

next_date <-
  tail(train$DATE, -1)

valid_transition <-
  as.integer(
    next_date - current_date
  ) == 1 &
  !is.na(current_state) &
  !is.na(next_state)

transition_counts <- table(
  Current = factor(
    current_state[valid_transition],
    levels = states
  ),
  Next = factor(
    next_state[valid_transition],
    levels = states
  )
)

transition_counts

# ------------------------------------------------------------
# 7. Transition probability matrix
# ------------------------------------------------------------

transition_matrix <-
  prop.table(
    transition_counts,
    margin = 1
  )

round(
  transition_matrix,
  3
)

rowSums(
  transition_matrix
)

# ------------------------------------------------------------
# 8. One-day-ahead prediction
# ------------------------------------------------------------

test_current <-
  head(test$State, -1)

test_actual <-
  tail(test$State, -1)

test_current_date <-
  head(test$DATE, -1)

test_next_date <-
  tail(test$DATE, -1)

valid_test <-
  as.integer(
    test_next_date - test_current_date
  ) == 1 &
  !is.na(test_current) &
  !is.na(test_actual)

test_current <-
  test_current[valid_test]

test_actual <-
  test_actual[valid_test]

predicted <- sapply(
  as.character(test_current),
  function(current) {

    probabilities <-
      transition_matrix[current, ]

    names(probabilities)[
      which.max(probabilities)
    ]
  }
)

# ------------------------------------------------------------
# 9. Markov model accuracy
# ------------------------------------------------------------

markov_accuracy <-
  mean(
    predicted ==
    as.character(test_actual)
  )

cat(
  "Markov Model Accuracy =",
  round(
    markov_accuracy * 100,
    2
  ),
  "%\n"
)

# ------------------------------------------------------------
# 10. Confusion matrix
# ------------------------------------------------------------

confusion_matrix <- table(
  Actual = factor(
    test_actual,
    levels = states
  ),
  Predicted = factor(
    predicted,
    levels = states
  )
)

confusion_matrix

# ------------------------------------------------------------
# 11. Persistence baseline
# ------------------------------------------------------------

persistence_prediction <-
  as.character(test_current)

persistence_accuracy <-
  mean(
    persistence_prediction ==
    as.character(test_actual)
  )

cat(
  "Persistence Baseline Accuracy =",
  round(
    persistence_accuracy * 100,
    2
  ),
  "%\n"
)

# ------------------------------------------------------------
# 12. Most-common-state baseline
# ------------------------------------------------------------

most_common_state <-
  names(
    which.max(
      table(train$State)
    )
  )

most_common_prediction <-
  rep(
    most_common_state,
    length(test_actual)
  )

most_common_accuracy <-
  mean(
    most_common_prediction ==
    as.character(test_actual)
  )

cat(
  "Most Common State =",
  most_common_state,
  "\n"
)

cat(
  "Most-Common-State Baseline Accuracy =",
  round(
    most_common_accuracy * 100,
    2
  ),
  "%\n"
)

# ------------------------------------------------------------
# 13. Compare prediction methods
# ------------------------------------------------------------

comparison <- data.frame(
  Method = c(
    "First-order Markov Chain",
    "Persistence Baseline",
    "Most-Common-State Baseline"
  ),
  Accuracy = c(
    markov_accuracy,
    persistence_accuracy,
    most_common_accuracy
  )
)

comparison$Accuracy_Percent <-
  round(
    comparison$Accuracy * 100,
    2
  )

print(comparison)
