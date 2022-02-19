#' Convert winning scores and labels into discoveries
#'
#' \code{fdp_sd} takes the output of \code{mirandom} and additional
#' statistical parameters to return the indices and winning scores of
#' hypotheses that were rejected.
#'
#' @param scores_and_labels An m x 2 matrix obtained via \code{mirandom}.
#' @param alpha An FDP threshold.
#' @param conf To control the FDP with \code{1 - conf} confidence.
#' @param c Determines the ranks of the target score that are considered
#' winning. Defaults to \code{c = 0.5} for single-decoy FDP-SD.
#' @param lambda Determines the ranks of the target score that are
#' considered losing. Defaults to \code{lambda = 0.5} for single-decoy FDP-SD.
#' @param procedure Takes a value of "standard" (for non-randomised FDP-SD) or
#' "coinflip" (for randomised FDP-SD).
#'
#' @return A list of 2 objects: the winning scores (\code{discoveries}) and
#' indices (\code{discoveries_ind}) of rejected hypotheses.
#' @export
#'
#' @examples
#' set.seed(123)
#' target_scores <- rnorm(200, mean = 1.5)
#' decoy_scores <- matrix(rnorm(600, mean = 0), ncol = 3)
#' scores <- cbind(target_scores, decoy_scores)
#' scores_and_labels <- mirandom(scores)
#' fdp_sd(scores_and_labels, alpha = 0.1, conf = 0.1)
fdp_sd = function(scores_and_labels,
                  alpha, conf,
                  c = 0.5, lambda = 0.5,
                  procedure = "standard") {

  if (c > lambda) {
    stop("Please choose c <= lambda.")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("Please choose alpha in (0, 1).")
  }
  if (conf <= 0 || conf >= 1) {
    stop("Please choose conf in (0, 1).")
  }
  if (procedure != "standard" && procedure != "coinflip") {
    stop("Please choose procedure to be [standard] or [coinflip].")
  }
  if (dim(scores_and_labels)[2] != 2) {
    stop("More than 2 columns detected in [scores_and_labels]. Are you using
       the correct matrix obtained from mirandom()?")
  }

  # Extract scores and labels
  W = scores_and_labels[, 1]
  L = scores_and_labels[, 2]
  n = length(W)
  n_count = sum(L != 0)
  counted_hyp = which(L != 0)

  # Permutation to break ties in W at random
  permW = sample(1:n_count, n_count, replace = FALSE)
  W_order = order(W[counted_hyp][permW], decreasing = TRUE)

  # First permute the winning scores and then order
  W_sort = W[counted_hyp][permW][W_order]
  L_sort = L[counted_hyp][permW][W_order]

  # Compute R(c, lambda) and i_0
  R_c_lam = (1 - lambda) / (1 - lambda + c)
  m_conf = ceiling(-log(1 / conf, base = 1 - R_c_lam) - 1e-13)
  i_0 = max(1, ceiling((m_conf - 1) / alpha - 1e-13))

  # Compute FDP-SD threshold
  if (i_0 > n_count) {
    threshold = 0
  } else {
    if (procedure == "coinflip") {
      delta_old = -1
      delta_new = -1
      bardelta_old = 0
      bardelta_new = 0
      omega_old = 1
      omega_new = 1
    }
    for (j in i_0:n_count) {
      d_j = sum((L_sort == -1)[1:j])
      for (d in 0:j) {
        k = floor((j - d) * alpha + 1e-13) + 1
        if (d != j & stats::pbinom(d, k + d, R_c_lam) - 1e-12 <= conf) {
          next
        } else if (
            d == j & stats::pbinom(d, k + d, R_c_lam) - 1e-12 <= conf
          ) {
          delta_new = j
          if (procedure == "coinflip") {
            p_floor = stats::pbinom(j, 1 + j, R_c_lam) - 1e-12
            p_ceil = 1
            omega_new = (p_ceil - conf) / (p_ceil - p_floor)
            if (bardelta_old == delta_new + 1) {
              bardelta_new = delta_new + 1
            } else if (delta_new > delta_old) {
              bardelta_new = delta_new + (stats::runif(1) > omega_new)
            } else {
              bardelta_new = delta_new +
                               (stats::runif(1) > omega_new / omega_old)
            }
          }
        } else {
          delta_new = d - 1
          if (procedure == "coinflip") {
            k = floor((j - d + 1) * alpha + 1e-13) + 1
            p_floor = stats::pbinom(d - 1, k + d - 1, R_c_lam) - 1e-12
            k = floor((j - d) * alpha + 1e-13) + 1
            p_ceil = stats::pbinom(d, k + d, R_c_lam)
            omega_new = (p_ceil - conf) / (p_ceil - p_floor)
            if (bardelta_old == delta_new + 1) {
              bardelta_new = delta_new + 1
            } else if (delta_new > delta_old) {
              bardelta_new = delta_new + (stats::runif(1) > omega_new)
            } else {
              bardelta_new = delta_new +
                               (stats::runif(1) > omega_new / omega_old)
            }
          }
        }
        if (procedure == "coinflip") {
          omega_old = omega_new
          delta_old = delta_new
          bardelta_old = bardelta_new
        }
        break
      }

      if (
        j == i_0 &
        d_j > ifelse(procedure == "coinflip", bardelta_new, delta_new)
      ) {
        threshold = 0
        break
      } else if (
        j != n_count &
        d_j <= ifelse(procedure == "coinflip", bardelta_new, delta_new)
      ) {
        next
      } else if (
        j == n_count &
        d_j <= ifelse(procedure == "coinflip", bardelta_new, delta_new)
      ) {
        threshold = j
        break
      } else {
        threshold = j - 1
        break
      }
    }
  }

  if(threshold > 0) {
    discoveries = which(L_sort[1:threshold] == 1)
    # Original indices of rejected hypotheses
    discoveries_ind = c(1:n)[counted_hyp][permW[W_order[discoveries]]]
    discoveries = W[discoveries_ind]
  } else {
    discoveries_ind = numeric(0)
    discoveries = numeric(0)
  }

  return(list(discoveries = discoveries, discoveries_ind = discoveries_ind))
}
