data_common$rank_difference <- data_common$score_cancer - data_common$score_healthy

# Calculate the sign of rank differences
data_common$sign_diff <- sign(data_common$rank_difference)

# Count the number of positive and negative differences
counts <- table(data_common$sign_diff)
num_positive <- counts["1"]
num_negative <- counts["-1"]

# Perform binomial test
binom_test <- binom.test(
  x = num_positive,
  n = num_positive + num_negative,
  p = 0.5,
  alternative = "two.sided"
)

# View the result
binom_test
