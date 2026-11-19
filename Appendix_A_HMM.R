\begin{lstlisting}[language=R]

# ============================================================
# HIDDEN MARKOV MODEL FOR DAILY RAINFALL
# Real Data: Wilmington International Airport
# NOAA Daily Weather Data, 2016-2025
# ============================================================

# Install these only once if needed
# install.packages("dplyr")
# install.packages("ggplot2")

library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Read and prepare Wilmington daily weather data
# ------------------------------------------------------------

wilmington <- read.csv(file.choose())

wilmington <- wilmington %>%
  mutate(
    date = as.Date(DATE),
    precipitation = PRCP,
    mean_temp = (TMAX + TMIN) / 2,
    mean_wind = AWND
  ) %>%
  arrange(date)

head(wilmington)

# ------------------------------------------------------------
# 2. Define rain occurrence
# ------------------------------------------------------------

wilmington$rain <-
  ifelse(wilmington$precipitation > 0, 1, 0)

head(wilmington)

table(wilmington$rain)

# ------------------------------------------------------------
# 3. Initial state probabilities
# ------------------------------------------------------------

pi <- c(0.70, 0.30)

# ------------------------------------------------------------
# 4. Transition probability matrix
# ------------------------------------------------------------

A <- matrix(
  c(0.85, 0.15,
    0.30, 0.70),
  nrow = 2,
  byrow = TRUE
)

rownames(A) <- c("Dry regime", "Wet regime")
colnames(A) <- c("Dry regime", "Wet regime")

A

# ------------------------------------------------------------
# 5. Emission probabilities
# ------------------------------------------------------------

rain_prob <- c(0.10, 0.80)

# ------------------------------------------------------------
# 6. Calculate emission probabilities
# ------------------------------------------------------------

emission_probability <- function(y, state) {

  p <- rain_prob[state]

  if (y == 1) {
    return(p)
  } else {
    return(1 - p)
  }
}

emission_probability(1, 1)
emission_probability(1, 2)

# ------------------------------------------------------------
# 7. Forward filtering algorithm
# ------------------------------------------------------------

forward_algorithm <- function(y, pi, A, rain_prob) {

  T <- length(y)
  K <- length(pi)

  alpha <- matrix(0, nrow = T, ncol = K)

  # Initial step
  for (j in 1:K) {

    emission <-
      ifelse(
        y[1] == 1,
        rain_prob[j],
        1 - rain_prob[j]
      )

    alpha[1, j] <- pi[j] * emission
  }

  # Normalize
  alpha[1, ] <- alpha[1, ] / sum(alpha[1, ])

  # Recursion
  for (t in 2:T) {

    for (j in 1:K) {

      emission <-
        ifelse(
          y[t] == 1,
          rain_prob[j],
          1 - rain_prob[j]
        )

      alpha[t, j] <-
        sum(alpha[t - 1, ] * A[, j]) * emission
    }

    alpha[t, ] <- alpha[t, ] / sum(alpha[t, ])
  }

  return(alpha)
}

# ------------------------------------------------------------
# 8. Apply the forward algorithm
# ------------------------------------------------------------

Y <- wilmington$rain

alpha <- forward_algorithm(
  Y,
  pi,
  A,
  rain_prob
)

colnames(alpha) <- c("Dry regime", "Wet regime")

head(alpha)

# ------------------------------------------------------------
# 9. Assign a hidden state to each day
# ------------------------------------------------------------

wilmington$state <-
  apply(alpha, 1, which.max)

table(wilmington$state)

# ------------------------------------------------------------
# 10. Interpret the hidden states
# ------------------------------------------------------------

state_summary <- wilmington %>%
  group_by(state) %>%
  summarise(
    Number_days = n(),
    Rain_frequency = mean(rain),
    Mean_precipitation = mean(precipitation),
    Mean_temperature = mean(mean_temp),
    Mean_wind = mean(mean_wind)
  )

state_summary

# ------------------------------------------------------------
# 11. Plot the inferred hidden regimes
# ------------------------------------------------------------

ggplot(
  wilmington,
  aes(
    x = date,
    y = precipitation,
    shape = factor(state)
  )
) +
  geom_point(size = 1.5) +
  labs(
    title =
      "Hidden Weather Regimes at Wilmington International Airport",
    subtitle =
      "Hidden Markov Model Classification",
    x = "Date",
    y = "Daily Precipitation",
    shape = "Hidden State"
  ) +
  theme_minimal()

# ------------------------------------------------------------
# 12. One-day-ahead hidden-state forecast
# ------------------------------------------------------------

current_prob <- alpha[nrow(alpha), ]

current_prob

next_state_prob <-
  as.numeric(current_prob %*% A)

names(next_state_prob) <-
  c("Dry regime", "Wet regime")

next_state_prob

# ------------------------------------------------------------
# 13. Forecast probability of rain tomorrow
# ------------------------------------------------------------

rain_tomorrow <-
  sum(next_state_prob * rain_prob)

cat(
  "Forecast probability of rain tomorrow =",
  round(rain_tomorrow, 3),
  "\n"
)
