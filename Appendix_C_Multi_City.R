# ============================================================
# ADDITIONAL LOCATIONS AND MULTI-LOCATION ANALYSIS
# Denver, Los Angeles, Seattle, and Miami
# NOAA Daily Weather Data, 2016-2025
# ============================================================

library(dplyr)

# ------------------------------------------------------------
# 1. Read the four additional datasets
# ------------------------------------------------------------

# Select dataset2.csv - Denver
denver <- read.csv(file.choose())

# Select dataset3.csv - Los Angeles
los_angeles <- read.csv(file.choose())

# Select dataset4.csv - Seattle
seattle <- read.csv(file.choose())

# Select dataset5.csv - Miami
miami <- read.csv(file.choose())

# Convert DATE to Date format
denver$DATE <- as.Date(denver$DATE)
los_angeles$DATE <- as.Date(los_angeles$DATE)
seattle$DATE <- as.Date(seattle$DATE)
miami$DATE <- as.Date(miami$DATE)

# Check station names
unique(denver$NAME)
unique(los_angeles$NAME)
unique(seattle$NAME)
unique(miami$NAME)

# Check dataset dimensions
dim(denver)
dim(los_angeles)
dim(seattle)
dim(miami)

# Check date ranges
range(denver$DATE)
range(los_angeles$DATE)
range(seattle$DATE)
range(miami$DATE)

# Check missing precipitation values
sum(is.na(denver$PRCP))
sum(is.na(los_angeles$PRCP))
sum(is.na(seattle$PRCP))
sum(is.na(miami$PRCP))


# ============================================================
# 2. CREATE FOUR WEATHER STATES
# ============================================================

# ------------------------------------------------------------
# 2.1 Create year variable
# ------------------------------------------------------------

denver$year <- as.integer(format(denver$DATE, "%Y"))
los_angeles$year <- as.integer(format(los_angeles$DATE, "%Y"))
seattle$year <- as.integer(format(seattle$DATE, "%Y"))
miami$year <- as.integer(format(miami$DATE, "%Y"))

# ------------------------------------------------------------
# 2.2 Split into training and testing data
# Training: 2016-2023
# Testing: 2024-2025
# ------------------------------------------------------------

denver_train <- denver %>%
  filter(year <= 2023)

denver_test <- denver %>%
  filter(year >= 2024)

la_train <- los_angeles %>%
  filter(year <= 2023)

la_test <- los_angeles %>%
  filter(year >= 2024)

seattle_train <- seattle %>%
  filter(year <= 2023)

seattle_test <- seattle %>%
  filter(year >= 2024)

miami_train <- miami %>%
  filter(year <= 2023)

miami_test <- miami %>%
  filter(year >= 2024)

# ------------------------------------------------------------
# 2.3 Heavy Rain thresholds
# 80th percentile of positive precipitation in training data
# ------------------------------------------------------------

denver_threshold <- as.numeric(
  quantile(
    denver_train$PRCP[denver_train$PRCP > 0],
    0.80,
    na.rm = TRUE
  )
)

la_threshold <- as.numeric(
  quantile(
    la_train$PRCP[la_train$PRCP > 0],
    0.80,
    na.rm = TRUE
  )
)

seattle_threshold <- as.numeric(
  quantile(
    seattle_train$PRCP[seattle_train$PRCP > 0],
    0.80,
    na.rm = TRUE
  )
)

miami_threshold <- as.numeric(
  quantile(
    miami_train$PRCP[miami_train$PRCP > 0],
    0.80,
    na.rm = TRUE
  )
)

cat(
  "Denver Heavy Rain threshold =",
  denver_threshold,
  "inches\n"
)

cat(
  "Los Angeles Heavy Rain threshold =",
  la_threshold,
  "inches\n"
)

cat(
  "Seattle Heavy Rain threshold =",
  seattle_threshold,
  "inches\n"
)

cat(
  "Miami Heavy Rain threshold =",
  miami_threshold,
  "inches\n"
)

# ------------------------------------------------------------
# 2.4 Define four weather states
# ------------------------------------------------------------

states <- c(
  "Dry",
  "Rain",
  "Heavy Rain",
  "Thunder-Reported Day"
)

classify_weather <- function(prcp, thunder, threshold) {

  state <- rep(
    NA_character_,
    length(prcp)
  )

  # Identify Thunder-Reported Day first
  storm_day <-
    !is.na(thunder) &
    thunder == 1

  state[storm_day] <- "Thunder-Reported Day"

  # Non-thunder-reported days
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

# ------------------------------------------------------------
# 2.5 Apply classification
# ------------------------------------------------------------

denver_train$State <- classify_weather(
  denver_train$PRCP,
  denver_train$WT03,
  denver_threshold
)

denver_test$State <- classify_weather(
  denver_test$PRCP,
  denver_test$WT03,
  denver_threshold
)

la_train$State <- classify_weather(
  la_train$PRCP,
  la_train$WT03,
  la_threshold
)

la_test$State <- classify_weather(
  la_test$PRCP,
  la_test$WT03,
  la_threshold
)

seattle_train$State <- classify_weather(
  seattle_train$PRCP,
  seattle_train$WT03,
  seattle_threshold
)

seattle_test$State <- classify_weather(
  seattle_test$PRCP,
  seattle_test$WT03,
  seattle_threshold
)

miami_train$State <- classify_weather(
  miami_train$PRCP,
  miami_train$WT03,
  miami_threshold
)

miami_test$State <- classify_weather(
  miami_test$PRCP,
  miami_test$WT03,
  miami_threshold
)

# ------------------------------------------------------------
# 2.6 Show state counts
# ------------------------------------------------------------

cat("\nDENVER TRAINING STATES\n")
print(
  table(
    denver_train$State,
    useNA = "ifany"
  )
)

cat("\nLOS ANGELES TRAINING STATES\n")
print(
  table(
    la_train$State,
    useNA = "ifany"
  )
)

cat("\nSEATTLE TRAINING STATES\n")
print(
  table(
    seattle_train$State,
    useNA = "ifany"
  )
)

cat("\nMIAMI TRAINING STATES\n")
print(
  table(
    miami_train$State,
    useNA = "ifany"
  )
)

cat("\nDENVER TESTING STATES\n")
print(
  table(
    denver_test$State,
    useNA = "ifany"
  )
)

cat("\nLOS ANGELES TESTING STATES\n")
print(
  table(
    la_test$State,
    useNA = "ifany"
  )
)

cat("\nSEATTLE TESTING STATES\n")
print(
  table(
    seattle_test$State,
    useNA = "ifany"
  )
)

cat("\nMIAMI TESTING STATES\n")
print(
  table(
    miami_test$State,
    useNA = "ifany"
  )
)


# ============================================================
# 3. TRANSITION MATRICES
# ============================================================

make_transition_matrix <- function(data) {

  current_state <-
    head(data$State, -1)

  next_state <-
    tail(data$State, -1)

  current_date <-
    head(data$DATE, -1)

  next_date <-
    tail(data$DATE, -1)

  valid <-
    as.integer(
      next_date - current_date
    ) == 1 &
    !is.na(current_state) &
    !is.na(next_state)

  counts <- table(
    Current = factor(
      current_state[valid],
      levels = states
    ),
    Next = factor(
      next_state[valid],
      levels = states
    )
  )

  probabilities <- prop.table(
    counts,
    margin = 1
  )

  list(
    counts = counts,
    probabilities = probabilities
  )
}

# ------------------------------------------------------------
# 3.1 Denver
# ------------------------------------------------------------

denver_transition <-
  make_transition_matrix(
    denver_train
  )

cat(
  "\nDENVER TRANSITION COUNT MATRIX\n"
)

print(
  denver_transition$counts
)

cat(
  "\nDENVER TRANSITION PROBABILITY MATRIX\n"
)

print(
  round(
    denver_transition$probabilities,
    3
  )
)

# ------------------------------------------------------------
# 3.2 Los Angeles
# ------------------------------------------------------------

la_transition <-
  make_transition_matrix(
    la_train
  )

cat(
  "\nLOS ANGELES TRANSITION COUNT MATRIX\n"
)

print(
  la_transition$counts
)

cat(
  "\nLOS ANGELES TRANSITION PROBABILITY MATRIX\n"
)

print(
  round(
    la_transition$probabilities,
    3
  )
)

# ------------------------------------------------------------
# 3.3 Seattle
# ------------------------------------------------------------

seattle_transition <-
  make_transition_matrix(
    seattle_train
  )

cat(
  "\nSEATTLE TRANSITION COUNT MATRIX\n"
)

print(
  seattle_transition$counts
)

cat(
  "\nSEATTLE TRANSITION PROBABILITY MATRIX\n"
)

print(
  round(
    seattle_transition$probabilities,
    3
  )
)

# ------------------------------------------------------------
# 3.4 Miami
# ------------------------------------------------------------

miami_transition <-
  make_transition_matrix(
    miami_train
  )

cat(
  "\nMIAMI TRANSITION COUNT MATRIX\n"
)

print(
  miami_transition$counts
)

cat(
  "\nMIAMI TRANSITION PROBABILITY MATRIX\n"
)

print(
  round(
    miami_transition$probabilities,
    3
  )
)


# ============================================================
# 4. TEST PREDICTION PERFORMANCE
# ============================================================

evaluate_markov <- function(
  train_data,
  test_data,
  transition_matrix
) {

  current_state <-
    head(test_data$State, -1)

  actual_state <-
    tail(test_data$State, -1)

  current_date <-
    head(test_data$DATE, -1)

  next_date <-
    tail(test_data$DATE, -1)

  valid <-
    as.integer(
      next_date - current_date
    ) == 1 &
    !is.na(current_state) &
    !is.na(actual_state)

  current_state <-
    current_state[valid]

  actual_state <-
    actual_state[valid]

  # Markov prediction
  predicted <- sapply(
    as.character(current_state),
    function(current) {

      probabilities <-
        transition_matrix[current, ]

      names(probabilities)[
        which.max(probabilities)
      ]
    }
  )

  # Markov accuracy
  markov_accuracy <- mean(
    predicted ==
    as.character(actual_state)
  )

  # Confusion matrix
  confusion <- table(
    Actual = factor(
      actual_state,
      levels = states
    ),
    Predicted = factor(
      predicted,
      levels = states
    )
  )

  # Persistence baseline
  persistence_prediction <-
    as.character(current_state)

  persistence_accuracy <- mean(
    persistence_prediction ==
    as.character(actual_state)
  )

  # Most-common-state baseline
  # The most common state is determined from training data
  most_common_state <- names(
    which.max(
      table(train_data$State)
    )
  )

  most_common_accuracy <- mean(
    most_common_state ==
    as.character(actual_state)
  )

  list(
    markov_accuracy = markov_accuracy,
    confusion = confusion,
    persistence_accuracy = persistence_accuracy,
    most_common_state = most_common_state,
    most_common_accuracy = most_common_accuracy
  )
}

# ------------------------------------------------------------
# 4.1 Denver
# ------------------------------------------------------------

denver_result <- evaluate_markov(
  denver_train,
  denver_test,
  denver_transition$probabilities
)

cat("\nDENVER RESULTS\n")

cat(
  "Markov Accuracy =",
  round(
    denver_result$markov_accuracy * 100,
    2
  ),
  "%\n"
)

cat(
  "Persistence Baseline =",
  round(
    denver_result$persistence_accuracy * 100,
    2
  ),
  "%\n"
)

cat(
  "Most Common State =",
  denver_result$most_common_state,
  "\n"
)

cat(
  "Most-Common-State Baseline =",
  round(
    denver_result$most_common_accuracy * 100,
    2
  ),
  "%\n"
)

print(
  denver_result$confusion
)

# ------------------------------------------------------------
# 4.2 Los Angeles
# ------------------------------------------------------------

la_result <- evaluate_markov(
  la_train,
  la_test,
  la_transition$probabilities
)

cat("\nLOS ANGELES RESULTS\n")

cat(
  "Markov Accuracy =",
  round(
    la_result$markov_accuracy * 100,
    2
  ),
  "%\n"
)

cat(
  "Persistence Baseline =",
  round(
    la_result$persistence_accuracy * 100,
    2
  ),
  "%\n"
)

cat(
  "Most Common State =",
  la_result$most_common_state,
  "\n"
)

cat(
  "Most-Common-State Baseline =",
  round(
    la_result$most_common_accuracy * 100,
    2
  ),
  "%\n"
)

print(
  la_result$confusion
)

# ------------------------------------------------------------
# 4.3 Seattle
# ------------------------------------------------------------

seattle_result <- evaluate_markov(
  seattle_train,
  seattle_test,
  seattle_transition$probabilities
)

cat("\nSEATTLE RESULTS\n")

cat(
  "Markov Accuracy =",
  round(
    seattle_result$markov_accuracy * 100,
    2
  ),
  "%\n"
)

cat(
  "Persistence Baseline =",
  round(
    seattle_result$persistence_accuracy * 100,
    2
  ),
  "%\n"
)

cat(
  "Most Common State =",
  seattle_result$most_common_state,
  "\n"
)

cat(
  "Most-Common-State Baseline =",
  round(
    seattle_result$most_common_accuracy * 100,
    2
  ),
  "%\n"
)

print(
  seattle_result$confusion
)

# ------------------------------------------------------------
# 4.4 Miami
# ------------------------------------------------------------

miami_result <- evaluate_markov(
  miami_train,
  miami_test,
  miami_transition$probabilities
)

cat("\nMIAMI RESULTS\n")

cat(
  "Markov Accuracy =",
  round(
    miami_result$markov_accuracy * 100,
    2
  ),
  "%\n"
)

cat(
  "Persistence Baseline =",
  round(
    miami_result$persistence_accuracy * 100,
    2
  ),
  "%\n"
)

cat(
  "Most Common State =",
  miami_result$most_common_state,
  "\n"
)

cat(
  "Most-Common-State Baseline =",
  round(
    miami_result$most_common_accuracy * 100,
    2
  ),
  "%\n"
)

print(
  miami_result$confusion
)


# ============================================================
# 5. FIVE-LOCATION MODEL COMPARISON
# ============================================================

five_city_comparison <- data.frame(

  Location = c(
    "Wilmington",
    "Denver",
    "Los Angeles",
    "Seattle",
    "Miami"
  ),

  Markov_Accuracy = c(
    66.44,
    72.05,
    89.86,
    66.34,
    71.51
  ),

  Persistence_Baseline = c(
    60.27,
    68.36,
    87.40,
    63.45,
    68.22
  ),

  Most_Common_State_Baseline = c(
    67.12,
    73.15,
    89.59,
    59.72,
    58.77
  )
)

print(
  five_city_comparison
)
