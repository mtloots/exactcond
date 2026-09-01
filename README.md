# exactcond

Exact conditional reference distributions for tables and for linear models whose structure was
selected from the same data. Pure C back-end, shared byte for byte with the Python package of the
same name.

## The one idea

A nuisance parameter is removed by conditioning on a statistic sufficient for it. What survives
carries no nuisance, so a reference drawn from it is exact at every sample size, with no appeal to
an asymptotic approximation and no bootstrap.

## Two conditionings

**Both margins fixed.** In an r by c table the row and column totals are sufficient for the
marginal parameters, so conditioning on both leaves the Fisher-Yates law. It is sampled exactly,
and by an argument needing no appeal to Patefield's algorithm: the conditional law given both
margins is precisely the law of the table formed by pairing a fixed vector of row labels with a
uniformly random permutation of the column labels.

**Direction only.** Under a null in which the tested structure does not enter the data-generating
process, a statistic invariant to the regression coefficients and to the error scale depends on the
data only through the residual direction, which is uniform on the sphere of the residual space
whatever those nuisance parameters are. The reference stays exact when the statistic is maximised
over candidates, and so remains valid after a spatial weight matrix has been selected by fit.

## Why the dispersion statistic

For nominal rows against an ordinal column, `stat = "dispersion"` is the between-row dispersion of
the mean score. Its components are the row mean scores themselves, so where the column score is a
signed difference the statistic is assembled out of the very contrasts being reported, and a reader
can follow the arithmetic from the table to the test. A rank test on such a variable answers a
different question, and answers it on data that are mostly ties.

## Example

```r
library(exactcond)
tab <- matrix(c(3, 228, 72, 10, 548, 67, 7, 276, 131), nrow = 3, byrow = TRUE)
ec_rc_test(tab, v = c(-1, 0, 1), B = 200000, seed = 4207)
```

## Notes on honesty

`ec_rc_test()` reports the largest null draw alongside the p-value, and says so when the observed
statistic exceeded every draw: that is a resolution floor, not a p-value, and the two should not be
written the same way. Every p-value comes with a Wilson interval for the Monte Carlo error in the
number of draws actually taken.

`ec_rc_two()` samples the paired law directly rather than convolving two separately estimated null
distributions. The convolution needs the observed score to be reachable under no association; for a
strongly associated table it is not, at any number of draws, and the failure arrives precisely when
both tables are strongly associated, which is the case in which the correct answer is a large
p-value.

## Reproducibility

Every replicate draws from its own stream, seeded from the base seed and the replicate index, so
results do not depend on the thread count and are identical across the R and Python fronts. A
byte-parity harness checks that agreement to fifteen significant figures over every entry point; it
lives alongside the Python front rather than in this repository, because it needs both trees present
to run.
