#' Convert target/decoy scores into winning scores and labels
#'
#' \code{mirandom} takes a collection of target and decoy scores for m
#' hypotheses and returns a winning score and label attached to each.
#' The argument \code{scores} must be organised in an m x (d + 1) matrix,
#' where d is the number of decoy scores. The first column of \code{scores}
#' must contain the target scores.
#'
#' @param scores An m x (d + 1) matrix where m is the number of hypothesis and
#' d is the number of decoy scores for each hypothesis. The first column of
#' \code{scores} are target scores and subsequent columns are decoy scores.
#' @param c Determines the ranks of the target score that are considered
#' winning. Defaults to \code{c = 0.5} for single-decoy FDP-SD.
#' @param lambda Determines the ranks of the target score that are
#' considered losing. Defaults to \code{lambda = 0.5} for single-decoy FDP-SD.
#'
#' @return A m x 2 matrix where m is the number of hypotheses. The first column
#' contains the winning scores and the second column contains the corresponding
#' labels.
#' @export
#'
#' @examples
#' target_scores <- rnorm(200, mean = 1.5)
#' decoy_scores <- matrix(rnorm(600, mean = 0), ncol = 3)
#' scores <- cbind(target_scores, decoy_scores)
#' mirandom(scores)
mirandom = function(scores, c = 0.5, lambda = 0.5) {

  if (c > lambda) {
    stop("Please choose c <= lambda.")
  }
  if (!is.matrix(scores)) {
    stop("Please input a valid matrix of target and decoy scores. Make sure
       the first column of [scores] are target scores and subsequent columns
       decoy scores.")
  }
  if (!(dim(scores)[2] > 1)) {
    stop("No decoy scores detected. Make sure the first column of the [scores]
       matrix are target scores and subsequent columns decoy scores.")
  }
  # Number of decoys
  n_p = ncol(scores) - 1
  if (
    !isTRUE(all.equal((((n_p + 1) * c)) %% 1, 0)) ||
    !isTRUE(all.equal((((n_p + 1) * lambda)) %% 1, 0))
  ) {
    stop("Please select both c and lambda of the form k/(d+1) where d is
       the number of decoy scores for each hypothesis and k is an integer
       between 1 and d (inclusive).")
  }

  ## Initialisation ##

  # Ranks of scores for each hypothesis
  rank_scores = t(apply(scores, 1, rank, ties.method = "random"))

  # Target scores
  obs_score = scores[, 1]

  # Decoy scores
  decoy_score = scores[, -1]

  # Number of hypotheses
  n = length(obs_score)

  # Storage for winning scores
  W = rep(0, n)

  # Storage for labels
  L = rep(0, n)

  ## Labelling Hypotheses ##

  # Calculate the empirical p-value of target scores
  pvals_obs = (n_p + 2 - rank_scores[, 1]) / (n_p + 1)

  # Which hypotheses are target-winning?
  ori_select = rep(FALSE, n)
  ori_select[which(pvals_obs - 1e-12 <= c)] = TRUE

  # Label target wins with 1
  L[ori_select] = 1

  # Which hypotheses are decoy-winning?
  decoy_select = rep(FALSE, n)
  decoy_select[which(pvals_obs - 1e-12 > lambda)] = TRUE

  # Label decoy wins with -1
  L[decoy_select] = -1

  # Which hypotheses are uncounted?
  uncounted = rep(FALSE, n)
  uncounted[!ori_select & !decoy_select] = TRUE

  # Label uncounted hypotheses with 0
  L[uncounted] = 0

  ## Computing Winning Scores ##

  # Target scores are winning for target-winning hypotheses
  W[ori_select] = scores[ori_select, 1]

  # Assign winning scores of uncounted hypotheses (selected uniformly at random)
  if (sum(uncounted) > 0) {

    # Location of uncounted hypotheses
    uncounted_loc = which(uncounted)

    # The number of winning ranks
    n_decoys_in_draw = floor(c * (n_p + 1) + 1e-12)

    # The minimum rank of winning ranks (d_1 - i_c + 1)
    min_rank = n_p + 2 - n_decoys_in_draw

    # The maximum rank of winning ranks (d_1)
    max_rank = n_p + 1

    # Randomly chosen ranks for each uncounted hypothesis
    rank_choice = as.integer(stats::runif(
                               sum(uncounted),
                               min = min_rank,
                               max = max_rank
                             ))

    # Matrix of ranks of uncounted hypotheses
    rank_uncounted_scores = as.matrix(rank_scores[uncounted_loc, ])

    # Assign winning scores
    for(i in 1:sum(uncounted))
    {
      unc_choice = which.min(abs(rank_uncounted_scores[i, ]  - rank_choice[i]))
      W[uncounted_loc[i]] = scores[uncounted_loc[i], unc_choice]
    }
  }

  # Indices of decoy-winning hypotheses
  decoy_select_loc = which(decoy_select)

  # Mirandom mapping (given we saw a decoy win)
  if(sum(decoy_select > 0)) {

    # If we have more than one decoy, the mapping is non-trivial
    if(n_p > 1) {

      # Number of winning ranks
      n_decoys_in_draw = floor(c*(n_p + 1) + 1e-12)

      # Number of losing ranks
      n_obs_ranks = n_p + 1 - ceiling((n_p + 1) * lambda - 1e-12)

      # Mirandom initialisation
      max_mapped_decoy_rank = pracma::zeros(n_obs_ranks, 1)
      max_mapped_decoy_prob = max_mapped_decoy_rank
      current_decoy_rank = n_p
      current_decoy_coverage = 0
      uni_decoy_coverage = n_obs_ranks / n_decoys_in_draw

      for (i in 1 : n_obs_ranks)  {
        max_mapped_decoy_rank[i] = current_decoy_rank

        if (current_decoy_coverage + 1 > uni_decoy_coverage) {
          max_mapped_decoy_prob[i] = uni_decoy_coverage - current_decoy_coverage
          remainder_obs_coverage_prob = 1 - max_mapped_decoy_prob[i]
          current_decoy_coverage = uni_decoy_coverage

        } else {
          max_mapped_decoy_prob[i] = 1
          current_decoy_coverage = 1 + current_decoy_coverage
          remainder_obs_coverage_prob = 0
        }

        if (current_decoy_coverage >= uni_decoy_coverage - 1e-10) {
          current_decoy_rank =
            current_decoy_rank -
            floor(remainder_obs_coverage_prob / uni_decoy_coverage + 1e-12) - 1

          current_decoy_coverage =
            remainder_obs_coverage_prob -
            floor(remainder_obs_coverage_prob / uni_decoy_coverage + 1e-12) *
            uni_decoy_coverage
        }
      }

      # Apply mirandom to rank of target scores for decoy-winning hypotheses
      obs_ranks = rank_scores[decoy_select, 1]
      rands = pracma::rand(sum(decoy_select), 1)

      mapped_obs_ranks =
        max_mapped_decoy_rank[obs_ranks] -
        ceiling(
          (rands - max_mapped_decoy_prob[obs_ranks]) /
            uni_decoy_coverage - 1e-12
        ) *
        (rands > max_mapped_decoy_prob[obs_ranks])

      # Adjust ranking of decoys
      rank_decoy_scores = as.matrix(rank_scores[decoy_select_loc, -1])
      rank_obs_scores = rank_scores[decoy_select_loc, 1]

      for(i in 1:sum(decoy_select)) {
        rank_decoy_scores[i, ][
          rank_decoy_scores[i, ] > rank_obs_scores[i]
        ] =
          rank_decoy_scores[i, ][
            rank_decoy_scores[i, ] > rank_obs_scores[i]
          ] - 1
      }

      # Assign winning scores
      for(i in 1:sum(decoy_select)) {
        dec_choice = which.min(
                       abs(rank_decoy_scores[i, ]  - mapped_obs_ranks[i])
                     )
        W[decoy_select_loc[i]] = scores[decoy_select_loc[i], dec_choice + 1]
      }
    } else {
      # If we only have one decoy, the winning score is the decoy score
      W[decoy_select_loc] = scores[decoy_select_loc, 2]
    }
  }

  scores_and_labels = matrix(c(W, L), ncol = 2, byrow = FALSE)
  return(scores_and_labels)
}
