delta_bounds = function(n, c, lambda, qs, confs, procedure = "standard",
                        delta_i.list) {
  cat("Beginning computation --", format(Sys.time(), "%a %b %d %X %Y"), "\n")

  R_c_lam                      <- (1 - lambda) / (1 - lambda + c)
  nqs                          <- length(qs)
  delta_table                  <- as.list(rep(0, nqs))
  names(delta_table)           <- qs
  for (i in 1:nqs) {
    delta_table[[i]]           <- matrix(0, nrow = n, ncol = length(confs))
    colnames(delta_table[[i]]) <- confs
  }

  for (v in seq_along(confs))  {
    cat(
      "Current confidence:", confs[v],
      "--",
      format(Sys.time(), "%a %b %d %X %Y"), "\n"
    )
    conf   <- confs[v]
    m_conf <- ceiling(-log(1 / conf, base = 1 - R_c_lam) - 1e-13)

    if (procedure == "coinflip") {
      omega_i.list              <- as.list(rep(0, nqs))
      omega_prime_i.list        <- as.list(rep(0, nqs))
      for (i in 1:nqs) {
        omega_i.list[[i]]       <- rep(1, n)
        omega_prime_i.list[[i]] <- rep(0, n)
      }
    }

    for (i in 1:nqs) {
      delta_i <- rep(-1, n)
      i_0     <- max(1, ceiling((m_conf - 1) / qs[i] - 1e-13))

      # Return no discoveries if i_0 > n
      if (i_0 > n) {
        delta_table[[i]][, v] <- delta_i
        break
      }

      for (j in i_0:n) {
        for (d in 0:j) {
          k_d                      <- floor((j - d) * qs[i] + 1e-13) + 1
          if (d != j & stats::pbinom(d, k_d + d, R_c_lam) - 1e-12 <= conf) {
            next
          } else if (d == j & stats::pbinom(d, k_d + d, R_c_lam) - 1e-12 <= conf) {
            delta_i[j]             <- d
            if (procedure == "coinflip") {
              p_floor_i            <- stats::pbinom(d, 1 + d, R_c_lam) - 1e-12
              p_ceil_i             <- 1
              omega_i              <- (p_ceil_i - conf) / (p_ceil_i - p_floor_i)
              omega_i.list[[i]][j] <- omega_i
            }
          } else {
            delta_i[j]             <- d - 1
            if (procedure == "coinflip") {
              k_d                  <- floor((j - d + 1) * qs[i] + 1e-13) + 1
              p_floor_i            <- stats::pbinom(d - 1, k_d + d - 1, R_c_lam) -
                1e-12
              k_d                  <- floor((j - d) * qs[i] + 1e-13) + 1
              p_ceil_i             <- stats::pbinom(d, k_d + d, R_c_lam)
              omega_i              <- (p_ceil_i - conf) / (p_ceil_i - p_floor_i)
              omega_i.list[[i]][j] <- omega_i
            }
          }
          break
        }
        if (procedure == "coinflip") {
          if (j >= i_0) {
            if (delta_i.list[[i]][j] == delta_i.list[[i]][j - 1]) {
              omega_prime_i.list[[i]][j] <- omega_i.list[[i]][j] /
                omega_i.list[[i]][j - 1]
            }
          }
        }
      }

      delta_table[[i]][, v] <- delta_i
    }
  }
  return(delta_table)
}
