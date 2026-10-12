## Submission

This is the first submission of exactcond, at version 0.2.0. Version 0.1.0 was prepared and checked
on win-builder in September but never submitted; 0.2.0 adds two exported functions,
`ec_diffnet_test()` (a permutation differential-network test) and `ec_multi_test()` (a single
permutation test for a battery of categorical items), with tests for both.

## Test environments

* local macOS (arm64), R 4.5, `R CMD check --as-cran`: 0 errors, 0 warnings, 2 notes.
* win-builder R-devel: see the result recorded below once it returns.
* GitHub Actions: Linux (R devel, release, oldrel-1), macOS and Windows (R release).

## R CMD check results

Two notes locally.

* New submission. Expected and unavoidable for a first submission. The DESCRIPTION URLs previously
  returned 404 and no longer do; the repository now exists and both the homepage and the issues page
  resolve.
* "Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy" and
  "Skipping checking math rendering: package 'V8' unavailable". Both are missing local tools rather
  than properties of the package; win-builder reported the HTML manual check as OK.

## Notes for the reviewer

The package compiles a small pure-C back-end with OpenMP. Each replicate draws from its own stream,
seeded from the base seed and the replicate index, so results do not depend on the number of
threads; there is a test asserting exactly that.

The correctness claim is exactness of the conditional reference distribution, so the tests do not
merely check that the code runs. The central test compares the sampler against the conditional law
obtained by COMPLETE ENUMERATION of all tables with the given margins, atom by atom, and a second
test compares a two by two case against `fisher.test`.
