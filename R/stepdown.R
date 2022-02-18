fdp_sd = function(scores_and_labels,
                  c, lambda, alpha, conf,
                  procedure = "standard") {
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
