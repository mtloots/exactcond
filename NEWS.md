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
