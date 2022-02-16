mdc_sd = function(c, lambda, scores, rank_scores, conf, alpha_range,
                  delta_i.list, n_count, L_sort, custom_scores, winning_given,
                  W_order, permW, n, counted_hyp, W) {

  # Compute winning scores and labels
  mirandom(scores, rank_scores, c, lambda, include_uncounted = FALSE)

  # Storage for the scores and indices of discoveries
  nqs = length(alpha_range)
  Discoveries = as.list(rep(0, nqs))
  Discoveries_ind = as.list(rep(0, nqs))

  R_c_lam = (1 - lambda) / (1 - lambda + c)
  m_conf = ceiling(-log(1 / conf, base = 1 - R_c_lam) - 1e-13)

  alpha_index = 1

  for (alpha in alpha_range) {

    delta_i.values = delta_i.list[[as.character(alpha)]]
    i_0 = max(1, ceiling((m_conf - 1) / alpha_range[alpha_index] - 1e-13))

    if (i_0 > n_count) {
      threshold = 0

    } else {

      for (j in i_0:n_count) {

        d_j = sum(L_sort[1:j])

        if (j == i_0 & d_j > delta_i.values[j]) {
          threshold = 0
          break

        } else if (j != n_count & d_j <= delta_i.values[j]) {
          next

        } else if (j == n_count & d_j <= delta_i.values[j]) {
          threshold = j
          break

        } else {
          threshold = j - 1
          break

        }
      }
    }

    if(threshold > 0) {

      if (!custom_scores || winning_given == FALSE) {
        Last_Disc = threshold
        Disc = which(L_sort[1:Last_Disc] == 0)

        # Find the original indices of the rejected hypotheses
        # before we did our reordering
        discoveries_ind = W_order[Disc]
        Discoveries_ind[[alpha_index]] = permW[discoveries_ind]
        Discoveries_ind[[alpha_index]] = c(1:n)[counted_hyp][
                                                 Discoveries_ind[[alpha_index]]
                                               ]

      } else {
        Last_Disc = threshold
        Discoveries_ind[[alpha_index]] = which(L_sort[1:Last_Disc] == 0)
      }

      Discoveries[[alpha_index]] = W[Discoveries_ind[[alpha_index]]]

    } else {
      Discoveries_ind[[alpha_index]] = numeric(0)
      Discoveries[[alpha_index]] = numeric(0)
    }

    alpha_index = alpha_index + 1
  }

  return(list(Discoveries = Discoveries, Discoveries_ind = Discoveries_ind))
}
