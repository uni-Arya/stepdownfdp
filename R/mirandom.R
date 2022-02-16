mirandom = function(scores, rank_scores, c, lambda, include_uncounted = TRUE,
                    custom_scores, winning_given)
{
  if (!custom_scores || winning_given == FALSE)
  {
    # Stores target and decoy scores
    obs_score = scores[1,]
    decoy_score = scores[-1,]

    # Number of decoys and hypotheses
    n = length(obs_score)
    n_p = nrow(scores) - 1

    # Stores the indices of the hypotheses which are: (a) target wins, (b) decoy wins
    ori_select = rep(FALSE, n)
    decoy_select = rep(FALSE, n)

    # Calculates the empirical p-value of the target scores
    pvals_obs = (n_p + 2 - rank_scores[1,]) / (n_p + 1)

    # Check which hypotheses are target, decoy wins and uncounted
    ori_select[which(pvals_obs - 1e-12 <= c)] = TRUE
    decoy_select[which(pvals_obs - 1e-12 > lambda)] = TRUE

    # Stores winning scores and labels for each hypothesis
    W = rep(0, n)
    L = rep(0, n)

    # If we have a target-winning hypothesis, the winning score is the target score
    W[ori_select] = scores[1, ori_select]

    # Stores indices of hypotheses in which a decoy won
    decoy_select_loc = which(decoy_select)

    # Label decoy wins with 1 (N.B., this is different to our theoretical work)
    L[decoy_select] = 1

    # Does our procedure care about uncounted scores?
    if (include_uncounted)
    {
      uncounted = rep(FALSE, n)
      uncounted[!ori_select & !decoy_select] = TRUE
      n_uncounted = sum(uncounted)

      # Assign winning scores of uncounted hypotheses (select, uniformly at random, a winning rank)
      if (n_uncounted > 0)
      {
        # Location of uncounted hypotheses
        uncounted_loc = which(uncounted)

        # The number of winning ranks
        n_decoys_in_draw = floor(c * (n_p + 1) + 1e-12)

        # The minimum rank of winning ranks (d_1 - i_c + 1)
        min_rank = n_p + 2 - n_decoys_in_draw

        # The maximum rank of winning ranks (d_1)
        max_rank = n_p + 1

        # Randomly chosen ranks for each uncounted hypothesis
        rank_choice = as.integer(stats::runif(n_uncounted, min = min_rank, max = max_rank))

        # Matrix of ranks of uncounted hypotheses
        rank_uncounted_scores = as.matrix(rank_scores[, uncounted_loc])

        # Assign winning scores
        for(i in 1:n_uncounted)
        {
          unc_choice = which.min(abs(rank_uncounted_scores[, i]  - rank_choice[i]))
          W[uncounted_loc[i]] = scores[unc_choice, uncounted_loc[i]]
        }
      }

      # Label uncounted wins with -1
      L[uncounted] = -1
    }

    # Mirandom mapping (given we saw a decoy win)
    if(sum(decoy_select > 0))
    {
      # If we have more than one decoy, the mapping is non-trivial
      if(n_p > 1)
      {
        # How does the mirandom mapping work?
        # Imagine that each winning rank is an empty bucket that can hold (1-lambda)/c liters of water,
        # and that each losing rank is a watering can that holds 1 liters of water.
        # Note that the total amount of water we have equals the total volume of all buckets.
        # We now play the following game: fill each bucket with the water from our watering cans, starting
        # from the first, non-empty bucket (the highest rank). If the bucket is filled with
        # water before the first watering can (the lowest rank) empties, move on to the next bucket
        # and continue to fill it with water. Otherwise, the bucket was not filled after we depleted
        # the watering can, so we use the next can, and so on, until all buckets are filled. The
        # question we wish to answer is: for each watering can, which bucket in our sequence did it attempt
        # to fill first, and what proportion of this bucket did it fill?

        n_decoys_in_draw = floor(c*(n_p + 1) + 1e-12)               # The number of empty buckets
        n_obs_ranks = n_p + 1 - ceiling((n_p + 1) * lambda - 1e-12) # The number of watering cans
        max_mapped_decoy_rank = pracma::zeros(n_obs_ranks, 1)               # The first part of our question
        max_mapped_decoy_prob = max_mapped_decoy_rank               # The second part of our question
        current_decoy_rank = n_p                                    # The bucket currently being filled (starting from the highest rank)
        current_decoy_coverage = 0                                  # How much of the bucket is filled so far?
        uni_decoy_coverage = n_obs_ranks / n_decoys_in_draw         # The volume of each bucket

        for (i in 1 : n_obs_ranks)  # For each watering can (starting from the lowest rank)
        {
          max_mapped_decoy_rank[i] = current_decoy_rank # Fill the current bucket

          if (current_decoy_coverage + 1 > uni_decoy_coverage) # If we pour all the water of this can into the bucket, the bucket will overfill
          {
            max_mapped_decoy_prob[i] = uni_decoy_coverage - current_decoy_coverage # So this is the amount it fills
            remainder_obs_coverage_prob = 1 - max_mapped_decoy_prob[i]             # And this is how much water is left in the can
            current_decoy_coverage = uni_decoy_coverage                            # The bucket is now full

          } else # Otherwise, we can pour all of the water of this can into the bucket
          {
            max_mapped_decoy_prob[i] = 1                        # And so it filled the bucket by 1L
            current_decoy_coverage = 1 + current_decoy_coverage # This is the amount of water in the bucket
            remainder_obs_coverage_prob = 0                     # And no water is left in the can
          }

          if (current_decoy_coverage >= uni_decoy_coverage - 1e-10) # Is the bucket full?
          {
            # The next watering can will fill which bucket now?
            current_decoy_rank = current_decoy_rank - floor(remainder_obs_coverage_prob / uni_decoy_coverage + 1e-12) - 1

            # And how much of that bucket was already filled by a previous watering can?
            current_decoy_coverage = remainder_obs_coverage_prob - floor(remainder_obs_coverage_prob / uni_decoy_coverage + 1e-12) * uni_decoy_coverage
          }
        }

        # Sanity check: did the final bucket overflow? Do we still have a bucket left to fill?
        if (n_decoys_in_draw > 0 && (current_decoy_coverage > 1e-10 || n_p - current_decoy_rank != n_decoys_in_draw))
        {
          stop('Tell Uri the mapping doesnt work...')
        }

        # Now, apply the mirandom map to every target score which lost in its corresponding hypothesis
        obs_ranks = rank_scores[1, decoy_select] # Ranks of the target score for each decoy-winning hypothesis
        rands = pracma::rand(sum(decoy_select), 1)       # Generate a vector of random numbers in [0, 1] to randomly assign the winning scores

        # Of all buckets that a particular watering can (partially) filled, pick one at random, with the probability of the choice being the
        # proportion of water from the can that went into each bucket
        mapped_obs_ranks = max_mapped_decoy_rank[obs_ranks] - ceiling((rands - max_mapped_decoy_prob[obs_ranks]) / uni_decoy_coverage - 1e-12) * (rands > max_mapped_decoy_prob[obs_ranks])

        # Note for the following section:
        # The losing target scores are mapped to: the (n_p-k)th rank among the decoy
        # scores for some 0 <= k <= n_p-1. Therefore, we need to first rank the decoy
        # scores appropriately (which is done in the following lines of code).

        rank_decoy_scores = as.matrix(rank_scores[-1,decoy_select_loc]) # Improper ranking (i.e., not a sequence 1, 2, 3, ..., r)
        rank_obs_scores = rank_scores[1,decoy_select_loc]

        # Make the ranking proper
        for(i in 1:sum(decoy_select))
        {
          rank_decoy_scores[,i][rank_decoy_scores[,i] > rank_obs_scores[i]] = rank_decoy_scores[,i][rank_decoy_scores[,i] > rank_obs_scores[i]] - 1
        }

        # Assigns the appropriate winning scores for decoy-winning hypotheses
        for(i in 1:sum(decoy_select))
        {
          dec_choice = which.min(abs(rank_decoy_scores[,i]  - mapped_obs_ranks[i]))
          W[decoy_select_loc[i]] = scores[dec_choice + 1, decoy_select_loc[i]] # +1 to account for observed score in first row
        }

      } else # If we use only one decoy, the winning score is the decoy score (the mirandom map is trivial)
      {
        W[decoy_select_loc] = scores[2, decoy_select_loc]
      }
    }


    if (include_uncounted)
    {
      # Use a permutation to randomly break ties in W
      permW = sample(1:n, n, replace = FALSE)
      W_order = order(W[permW], decreasing = TRUE)

      # Order the scores and corresponding labels in descending fashion
      W_sort = W[permW][W_order]
      L_sort = L[permW][W_order]

    } else
    {
      # Indices of counted hypotheses
      counted_hyp = as.logical(ori_select + decoy_select)

      # The number of counted hypotheses
      n_count = sum(counted_hyp)

      # Use a permutation to randomly break ties in W
      permW = sample(1:n_count, n_count, replace = FALSE)
      W_order = order(W[counted_hyp][permW], decreasing = TRUE)

      # Order the scores and corresponding labels in descending fashion
      W_sort = W[counted_hyp][permW][W_order]
      L_sort = L[counted_hyp][permW][W_order]
    }

  } else # Use custom scores
  {
    # Define winning scores and labels
    W = scores
    L_sort = W
    decoy_wins = which(W < 0)
    target_wins = which(W > 0)
    uncounted_wins = which(W == 0)
    L_sort[decoy_wins] = 1
    L_sort[target_wins] = 0

    # Are we including uncounted scores and labels?
    if (include_uncounted)
    {
      L_sort[uncounted_wins] = -1
    } else
    {
      W = W[!uncounted_wins]
      L_sort = L_sort[!uncounted_wins]
      n_count = length(decoy_wins) + length(target_wins)
    }

  }

  # Assign variables to parent frame
  assign("W", W, pos = parent.frame())
  assign("L_sort", L_sort, pos = parent.frame())
  if (!custom_scores || winning_given == FALSE)
  {
    assign("permW", permW, pos = parent.frame())
    assign("W_order", W_order, pos = parent.frame())
  }
  if (!include_uncounted)
  {
    assign("n_count", n_count, pos = parent.frame())
    assign("counted_hyp", counted_hyp, pos = parent.frame())
  }
}
