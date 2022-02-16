#### Parameter generation ####


# Generates parameters as described in mdc.R
# calibration = 0 => calibrated scores
# calibration = 1 => non-calibrated scores

# There's also calibration = 2 (experimental)
# 90% of mu_i are 0, while 10% of mu_i are 150
# Variance sigma_i^2 = 1 is fixed for all hypotheses, separation nu_i is like in non-calibration

Parameter_generation = function(n, k, separation, mu_mu, sigma_mu, omega, const_par_nu, calibration = 0)
{
  # Experiment
  if (calibration == 2)
  {
    nu = 1 + stats::rexp(k, separation)
    mu = rep(0, n)
    mu[sample(1:n, n, replace = FALSE)] = c(rep(0, n*0.9), rep(150, n*0.1))
    sigma = rep(1, n)

  } else if (calibration == 1) # Calibration
  {
    nu = rep(const_par_nu, k)
    mu = rep(0, n)
    sigma = rep(1, n)

  } else # Non-calibration
  {
    nu = 1 + stats::rexp(k, separation)
    mu = stats::rnorm(n, mu_mu, sigma_mu)
    sigma = 1 + stats::rexp(n, 1/omega)
  }

  # Result stored in a list where each item contains a vector (lengths: k, n, n)
  return(list(nu = nu, mu = mu, sigma = sigma))
}


#### Target score generation ####


# Generates the observed scores by sampling them from Normal(mu_i + nu_i, sigma_i^2)
# for the alternatives or from a Normal(mu_i, sigma_i^2) for a true null.
Target_score_generation = function(n, mu, sigma)
{
  scores = stats::rnorm(n, mu, sigma)

  # Result stored in a vector of length n
  return(scores)
}


#### Decoy score generation ####


# Generates n_p decoy scores for each hypothesis, sampling from a Normal(mu_i, sigma_i^2)
Decoy_Score_generation = function(n, n_p, mu, sigma)
{
  decoy_score = matrix(stats::rnorm(n*n_p, mu, sigma), ncol = n, nrow = n_p, byrow = TRUE)

  # Result stored in a n_p x n matrix (each column represents the decoy scores for a single hypothesis)
  return(decoy_score)
}


#### Generate target and decoy scores for simulation ####


generate_scores = function(custom_scores, seeds, i_run, static_alt_score,
                           n, k, separation, mu_mu, sigma_mu, omega, n_decoys,
                           seed, parallel_compute, n_cores, winning_given,
                           score_list, n_vec, sample_known, calibration,
                           correct_rejections_list, const_par_nu)
{
  if (!custom_scores)
  {
    # Fixed seed for each run
    set.seed(strtoi(seeds[[i_run]], 16L))

    # Generate parameters for null distributions
    # Only generate new parameters if we don't want target scores of false nulls to be fixed
    if(static_alt_score == FALSE || i_run == 1)
    {
      # Generate Parameters (mu_i, sigma_i, nu_i)
      parameters = Parameter_generation(n, k, separation, mu_mu, sigma_mu, omega, const_par_nu, calibration)
      mu = parameters$mu
      sigma = parameters$sigma
      nu = parameters$nu

      # Assigns means and variances of true (null) and false (alt) hypotheses
      null_mu = mu[1:(n-k)]
      alt_mu = mu[(n-k+1):n] + nu
      null_sigma = sigma[1:(n-k)]
      alt_sigma = sigma[(n-k+1):n]
    }

    # Randomly assigns the indices of false nulls
    falseNullInds[,i_run] = sample(1:n, k, replace = FALSE)

    # Indicates which hypotheses are truly null (TRUE => true null)
    trueNullInds = rep(TRUE,n)
    trueNullInds[falseNullInds[,i_run]] = FALSE

    # If we have no false nulls, all scores are generated from a N(mu_i, sigma_i^2)
    if(k == 0)
    {
      obs_scores = Target_score_generation(n, mu, sigma)

      # Store mu_i and sigma_i for decoy score generation
      imu[,i_run] = mu
      isigma[,i_run] = sigma

    } else # We have some false nulls
    {
      # If, for each run, we want the target scores of false nulls to be fixed, generate them only once
      if(static_alt_score == FALSE || i_run == 1)
      {
        # Generate target scores for false nulls
        alt_scores = Target_score_generation(k, alt_mu, alt_sigma)
      }

      # Generate target scores for true nulls
      null_scores = Target_score_generation(n-k, null_mu, null_sigma)

      # Create a vector of target scores
      obs_scores = rep(0,n)
      obs_scores[falseNullInds[,i_run]] = alt_scores
      obs_scores[-falseNullInds[,i_run]] = null_scores

      # Store mu_i for decoy score generation
      imu[,i_run][falseNullInds[,i_run]] = alt_mu - nu
      imu[,i_run][-falseNullInds[,i_run]] = null_mu

      # Store sigma_i for decoy score generation
      isigma[,i_run][falseNullInds[,i_run]] = alt_sigma
      isigma[,i_run][-falseNullInds[,i_run]] = null_sigma
    }

    # Generate (matrix of) decoy scores
    total_decoy_scores = Decoy_Score_generation(n, n_decoys, imu[,i_run], isigma[,i_run])

    # Can be removed
    set.seed(strtoi(seed, 16L))

    # Turn the collection of scores into a matrix where each column represents the scores for a single hypothesis,
    # the first of which being the target
    scores = rbind(obs_scores, total_decoy_scores)

    # Rank the scores of each column, breaking ties at random
    if (parallel_compute) {
      n_clusters = min(n_cores, 2)
      cl = parallel::makeCluster(n_clusters)
      rank_scores = parallel::parApply(cl, scores, 2, rank, ties.method = "random")
      parallel::stopCluster(cl)
    } else {
      rank_scores = apply(scores, 2, rank, ties.method = "random")
    }

    # Sanity check (we should never run into this problem)
    obs_scores_pvals = (n_decoys + 2 - rank_scores[1,]) / (n_decoys + 1)
    if(any(obs_scores_pvals > 1 | obs_scores_pvals < 0))
    {
      stop("Issue with p-values.")
    }

  } else if (winning_given == TRUE) # We are given winning scores/labels
  {
    # Load scores and their ranks
    scores = score_list[[i_run]]
    rank_scores = NA

    # Do we know which of our scores correspond to true nulls?
    if (sample_known == TRUE)
    {
      correct_rejections = correct_rejections_list[[i_run]]
    }
  } else if (winning_given == FALSE)
  {
    # Load scores and their ranks
    scores = score_list[[i_run]]

    if (parallel_compute) {
      n_clusters = min(n_cores, 2)
      cl = parallel::makeCluster(n_clusters)
      rank_scores = parallel::parApply(cl, scores, 2, rank, ties.method = "random")
      parallel::stopCluster(cl)
    } else {
      rank_scores = apply(scores, 2, rank, ties.method = "random")
    }

    # Do we know which of our scores correspond to true nulls?
    if (sample_known == TRUE)
    {
      correct_rejections = correct_rejections_list[[i_run]]
    }
  }

  if (!custom_scores) {
    return(list(scores = scores, rank_scores = rank_scores, trueNullInds = trueNullInds))
  }

  if (sample_known == TRUE) {
    return(list(scores = scores, rank_scores = rank_scores, n = n, correct_rejections = correct_rejections))
  } else {
    return(list(scores = scores, rank_scores = rank_scores, n = n))
  }

}
