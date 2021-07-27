library(tidyverse)

reps <- 1000

prior_df <- tibble(
  seed = sample(0:1e8, size = reps, replace = FALSE),
  sample_sizes = rep("12-20-100-110", reps),
  mu = rep(3e-6, reps),
  c = rep(1e-6, reps),
  loci = rep(2e7, reps),
  neg_mean = -runif(reps, 0, 0.05),
  neg_shape = runif(reps, 0, 1),
  pos_mean = runif(reps, 0, 0.05),
  pos_shape = runif(reps, 0, 1),
  neg_prop = runif(reps, 0.90, 1),
  pos_prop = 1-neg_prop,
  na = as.integer(runif(reps, 10, 2000)), #!
  nb = as.integer(runif(reps, 10, 2000)), #!
  n0 = as.integer(runif(reps, 55, 2000)), #!
  tb = as.integer(runif(reps, 0, na)),
  t0 = as.integer(runif(reps, 0, na)),
)
readr::write_csv(x = prior_df, path = "prior_df.csv")
