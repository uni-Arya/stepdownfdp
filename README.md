
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stepdownfdp

<!-- badges: start -->
<!-- badges: end -->

This package provides a step-down procedure for controlling the False
Discovery Proportion (FDP) in a competition-based setup. This includes
target-decoy competition (TDC) in computational mass spectrometry and
the knockoff construction in regression. FDP control (also referred to
as FDX) is probabilistic in nature: given a prespecified FDP tolerance
*α* ∈ (0,1) and a desired confidence level 1 − *γ* the procedure reports
a list of discoveries for which the FDP is  ≤ *α* with probability
 ≥ 1 − *γ*.

## Installation

You can install the development version of stepdownfdp from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("uni-Arya/stepdownfdp")
```

## Usage

### With a single decoy

For use in a single-decoy experiment, the user must input a collection
of target/decoy scores for each hypothesis, an FDP threshold
*α* ∈ (0,1), and a confidence level *γ* ∈ (0,1).

First use `mirandom()` to convert target/decoy scores into winning
scores and labels. The argument `scores` must be an *m* × 2 matrix,
where *m* is the number of hypotheses. Make sure that the first column
of `scores` corresponds to target scores.

``` r
library(stepdownfdp)
set.seed(123)
target_scores     <- rnorm(200, mean = 1.75)
decoy_scores      <- rnorm(200, mean = 0)
scores            <- cbind(target_scores, decoy_scores)
scores_and_labels <- mirandom(scores)
head(scores_and_labels)
#>          [,1] [,2]
#> [1,] 2.198810   -1
#> [2,] 1.519823    1
#> [3,] 3.308708    1
#> [4,] 1.820508    1
#> [5,] 1.879288    1
#> [6,] 3.465065    1
```

Pass the output of `mirandom()` into `fdp_sd()` to return the winning
scores and indices of hypotheses that were rejected. The output of
`fdp_sd()` is in decreasing order of winning scores.

``` r
results  <- fdp_sd(scores_and_labels, alpha = 0.1, conf = 0.1)
W_scores <- results$discoveries
indices  <- results$discoveries_ind
W_scores
#>  [1] 4.991040 3.937333 3.918956 3.878452 3.850109 3.800085 3.747213 3.659104
#>  [9] 3.593862 3.536913 3.465065 3.308708 3.282611 3.266471 3.194551 3.118602
#> [17] 3.110652 3.013185 3.003815 2.974082 2.957962 2.898808 2.881337 2.859920
#> [25] 2.846839 2.802711 2.775571 2.755739 2.743504 2.726973 2.672267 2.668997
#> [33] 2.645126 2.628133 2.587787 2.571581
indices
#>  [1] 164  97  44 174 149  70 196 139 125  16   6   3  98  56 131  54  95 182  30
#> [20]  11  45  90 136 187  87 161  76  73  91 159  69 110  33  34  27  35
```

### With multiple decoys

For use in a multiple-decoy experiment, the user must also input
parameters *c* ≤ *λ* of the form *k*/(*d*+1) where *d* is the number of
decoy scores used for each hypothesis and 1 ≤ *k* ≤ *d* is an integer.
As an example, if we compute 3 decoy scores for each hypothesis, we may
take *c* and *λ* to be 1/4, 1/2, or 3/4, subject to *c* ≤ *λ*. The value
of *c* determines the ranks in which the target score is considered
“winning.” E.g., if *c* = 1/4, *H*<sub>*i*</sub> is labelled as a target
win whenever its corresponding target score is the highest ranked score
among all targets/decoys for that hypothesis. Similarly, *λ* determines
when the target score is “losing.”

The argument `scores` of `mirandom()` must now be an *m* × (*d*+1)
matrix, where *m* is the number of hypotheses and *d* is the number of
decoy scores for each. As in the single-decoy case, the first column
must consist of target scores.

``` r
library(stepdownfdp)
set.seed(123)
target_scores     <- rnorm(200, mean = 1.75)
decoy_scores      <- matrix(rnorm(600, mean = 0), ncol = 3)
scores            <- cbind(target_scores, decoy_scores)
scores_and_labels <- mirandom(scores, c = 0.25, lambda = 0.75)
head(scores_and_labels)
#>          [,1] [,2]
#> [1,] 2.198810    0
#> [2,] 1.519823    1
#> [3,] 3.308708    1
#> [4,] 1.820508    1
#> [5,] 1.879288    1
#> [6,] 3.465065    1
```

Pass the output into `fdp_sd()`, making sure to specify the same *c* and
*λ* values used in `mirandom()`.

``` r
results  <- fdp_sd(scores_and_labels, alpha = 0.1, conf = 0.1, c = 0.25, lambda = 0.75)
W_scores <- results$discoveries
indices  <- results$discoveries_ind
W_scores
#>  [1] 4.991040 3.937333 3.918956 3.878452 3.850109 3.800085 3.747213 3.659104
#>  [9] 3.593862 3.536913 3.465065 3.308708 3.282611 3.266471 3.194551 3.118602
#> [17] 3.110652 3.013185 3.003815 2.974082 2.957962 2.898808 2.881337 2.859920
#> [25] 2.846839 2.802711 2.775571 2.755739 2.743504 2.726973 2.672267 2.668997
#> [33] 2.645126 2.628133 2.587787 2.571581 2.537739 2.529965 2.519042 2.504054
#> [41] 2.489948 2.451784 2.451356 2.438640 2.437917 2.394377 2.386570
indices
#>  [1] 164  97  44 174 149  70 196 139 125  16   6   3  98  56 131  54  95 182  30
#> [20]  11  45  90 136 187  87 161  76  73  91 159  69 110  33  34  27  35 151  49
#> [39] 152 189 138 141  19  36 148  84 167
```

### Randomised versions

The user can also invoke the randomised version of both single-decoy and
multiple-decoy procedures by passing `procedure = "coinflip"` into
`fdp_sd()`.

``` r
results  <- fdp_sd(scores_and_labels, alpha = 0.1, conf = 0.1, c = 0.25, lambda = 0.75, 
                   procedure = "coinflip")
W_scores <- results$discoveries
indices  <- results$discoveries_ind
W_scores
#>   [1] 4.9910399 3.9373330 3.9189560 3.8784519 3.8501089 3.8000847 3.7472134
#>   [8] 3.6591036 3.5938620 3.5369131 3.4650650 3.3087083 3.2826106 3.2664706
#>  [15] 3.1945509 3.1186023 3.1106524 3.0131852 3.0038149 2.9740818 2.9579620
#>  [22] 2.8988076 2.8813372 2.8599203 2.8468390 2.8027115 2.7755714 2.7557385
#>  [29] 2.7435039 2.7269734 2.6722675 2.6689966 2.6451257 2.6281335 2.5877870
#>  [36] 2.5715811 2.5377388 2.5299651 2.5190422 2.5040538 2.4899475 2.4517843
#>  [43] 2.4513559 2.4386403 2.4379168 2.3943765 2.3865697 2.3579643 2.3507088
#>  [50] 2.3346137 2.3129895 2.3039177 2.2983970 2.2694072 2.2668620 2.2478505
#>  [57] 2.2109162 2.2015041 2.1982098 2.1865235 2.1851815 2.1764642 2.1507715
#>  [64] 2.1352804 2.1296395 2.1189645 2.1098138 2.0822026 2.0817820 2.0604807
#>  [71] 2.0535286 2.0511534 2.0482276 2.0068837 2.0033185 1.9887317 1.9853866
#>  [78] 1.9659416 1.9644453 1.9313035 1.9033731 1.8792877 1.8738542 1.8676466
#>  [85] 1.8606827 1.8556762 1.8445835 1.8347373 1.8279608 1.8205084 1.8152930
#>  [92] 1.8030042 1.7912329 1.7877884 1.7557642 1.7214532 1.7159327 1.7071295
#>  [99] 1.7049723 1.6944380 1.6880883 1.6786919 1.6666309 1.6111086 1.5528241
#> [106] 1.5420827 1.5320251 1.5295134 1.5198225 1.5142996 1.5137204 1.5033081
#> [113] 1.4939078 1.4878025 1.4696047 1.4549285 1.4440373 1.4253141 1.4024574
#> [120] 1.4003496 1.3793400 1.3775612 1.3697735 1.3695290 1.3471152 1.3331424
#> [127] 1.3043380 1.2916347 1.2833446 1.2662194 1.2594426 1.2507080 1.2190935
#> [134] 1.1941589 1.1249607 1.1220939 1.1092940 1.0980501 1.0619914 1.0552930
#> [141] 1.0407992 1.0395934 1.0086639 0.9650955 0.8844871 0.8546366 0.7983814
#> [148] 0.7416234 0.7008230 0.5292823 0.4629695
indices
#>   [1] 164  97  44 174 149  70 196 139 125  16   6   3  98  56 131  54  95 182
#>  [19]  30  11  45  90 136 187  87 161  76  73  91 159  69 110  33  34  27  35
#>  [37] 151  49 152 189 138 141  19  36 148  84 167 112 197  58 157  37  92 115
#>  [55] 169  17   7 132  67 179  88  31  13  82  61 170  12 153  86 178  66 116
#>  [73] 166 102  51  93 127  60 191  79  28   5  59 121  14 117 193 188 128   4
#>  [91] 172  68 133 177  81  52 173  53 106 114  38 130  50  80 186  42  22  85
#> [109]   2  99 185 103 124 142 156  32  39 192 104 183  83 158 109  40  47 165
#> [127]  10 180  48 168 123 190 146  15  25  94 118 126  75  41  74 101 175 107
#> [145] 184 194 105 154 162  78 150
```
