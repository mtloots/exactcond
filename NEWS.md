# exactcond 0.2.0

* `ec_diffnet_test()`: a differential network without a Gaussian assumption, following the
  permutation route of NetCoMi's `diffnet` and the max-over-edges statistic of the Network
  Comparison Test, with exchangeability blocks in the sense of Winkler et al. (2014). The edge statistic is
  the heterogeneity of the Fisher-transformed correlations across groups, the same quantity the
  Gaussian theory uses, referred to the law obtained by permuting the group labels instead of to a
  chi-square; on ranks the procedure is distribution-free. One p-value for the whole network (sum or
  maximum over edges) with the per-edge p-values from the same permutations for a closed follow-up,
  so a p(p-1)/2-edge network pays no multiplicity correction, and a cluster permutation for
  interviewer or ward effects. Measured cost against the Gaussian route on Gaussian data: about five
  points of power. Measured gain off it: the Gaussian null rejects a true null 28 per cent of the
  time on t(3) data and 29 per cent under clustering, where this holds 5 and 4.
* `ec_multi_test()`: one permutation test for a battery of categorical items across groups. This is
  Pesarin's nonparametric combination with the direct (sum) and Tippett (maximum) combining
  functions, and the maximum with per-item follow-up is a step of Westfall and Young's maxT;
  the sum or the maximum over items of a per-item Pearson or score-dispersion statistic referred to
  the law obtained by permuting the group labels, across respondents or across clusters. One p-value
  for the whole battery, no multiplicity correction, and per-item p-values from the same
  permutations for a closed follow-up. Provided because no package in this programme's toolchain
  offered these with exchangeability blocks in one place; on the prevention-behaviour paper's data
  the maximum-statistic form matched or beat a Bonferroni-corrected Fisher battery at its own size.

# exactcond 0.1.0

First release.

* `ec_rc_test()` and `ec_rc_null()`: exact conditional inference in a two-way table given both
  margins, sampled by random permutation of the column labels. Statistics are Goodman's
  linear-by-linear association, the dispersion of the row mean scores for nominal rows against an
  ordinal column, Pearson's statistic and the likelihood ratio.
* `ec_rc_two()`: exact conditional test that two tables share a linear-by-linear parameter.
* `ec_sphere_maxquad()`: exact null for a quadratic form maximised over candidate structures,
  valid after selection by fit.
* `ec_pvalue()`: tail probabilities with a mid-p companion and a Wilson interval for the Monte
  Carlo error, so an estimated p-value is not quoted as though it were exact arithmetic.
