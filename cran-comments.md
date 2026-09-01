## Test environments

* local macOS 15 (arm64), R 4.5.1, `R CMD check --as-cran`: 0 errors, 0 warnings, 2 notes.

Nothing else has been run yet. In particular there has been no win-builder run, no r-devel run and
no GitHub Actions matrix; those must be done before submission and this file updated to say what
they actually reported, rather than what they would be expected to report.

## R CMD check results

Two notes, both of them properties of the checking environment rather than of the package.

* New submission, and the URLs in DESCRIPTION return 404 because the repository has not yet been
  created. To be resolved before submission, either by creating it or by removing the fields.
* "Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy" and
  "Skipping checking math rendering: package 'V8' unavailable". Both are missing local tools.

## Notes for the reviewer

The package compiles a small pure-C back-end with OpenMP. Each replicate draws from its own stream,
seeded from the base seed and the replicate index, so results do not depend on the number of
threads; there is a test asserting exactly that.

The correctness claim is exactness of the conditional reference distribution, so the tests do not
merely check that the code runs. The central test compares the sampler against the conditional law
obtained by COMPLETE ENUMERATION of all tables with the given margins, atom by atom, and a second
test compares a two by two case against `fisher.test`.
