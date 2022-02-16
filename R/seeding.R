#### Seeding ####

seeding = function(n_runs) {
  seed = "1a50a8f8" # If 0 the seed will be randomized based on the time, otherwise the hexadecimal string will be used

  # If seed = 0, randomly select a new seed
  if(is.numeric(seed) && seed == 0){ # Distinguishes between seed = "0" and seed = 0; the former case is when we want the seed to be precisely 0

    seed_dec = as.integer(stats::runif(1, 1, .Machine$integer.max))
    seed = as.hexmode(seed_dec)

  }

  # Either use the seed from [### Other parameters ###] or the randomly generated seed above
  set.seed(strtoi(seed,16L))

  # Save the RNG state for the generation of the RNG for each run (see below)
  oldseed = .Random.seed

  # For each run, we generate a new RNG state
  # Used for reproducibility
  seeds = as.list(rep(0, n_runs))

  for(i in 1:n_runs){

    seed_dec = as.integer(stats::runif(1, 1, .Machine$integer.max))
    seeds[[i]] = as.hexmode(seed_dec)

  }

  return(list(seed = seed, seeds = seeds))
}
