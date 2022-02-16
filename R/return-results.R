run_methods = function(test_methods, seeds, i_run, scores, rank_scores, fixed_c,
                       c_lam, conf, qs, delta_i.list, nqs, custom_scores, FDPs,
                       trueNullInds, nTDs, sample_known, correct_rejections,
                       num_discoveries)
{
  if('mdc_sd' %in% test_methods)
  {
    set.seed(strtoi(seeds[[i_run]], 16L))

    # Compute, for each alpha, the rejection set (as a list indexed by alpha)
    result = mdc_sd(scores, rank_scores, fixed_c, max(fixed_c, c_lam), conf, qs, delta_i.list)
    R_list = result$Discoveries_ind

    # Store important quantities
    for(iq in 1:nqs)
    {
      R = R_list[[iq]]

      if (!custom_scores)
      {
        FDPs$mdc_sd[iq, i_run] <<- sum(trueNullInds[R]) / max(length(R), 1)
        nTDs$mdc_sd[iq, i_run] <<- sum(trueNullInds[R] == FALSE)
      } else if (sample_known)
      {
        FDPs$mdc_sd[iq, i_run] <<- sum(correct_rejections[R] == FALSE) / max(length(R), 1)
        nTDs$mdc_sd[iq, i_run] <<- sum(correct_rejections[R] == TRUE) / sum(correct_rejections)
      }
      num_discoveries$mdc_sd[iq, i_run] <<- pracma::numel(R)
    }

  }
}
