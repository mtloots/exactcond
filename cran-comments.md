## Test environments

* local macOS 15 (arm64), R 4.5.1, `R CMD check --as-cran`: 0 errors, 0 warnings, 2 notes.
* win-builder R-devel and R-release: submitted 1 September 2026, RESULTS NOT YET RECEIVED.

No GitHub Actions matrix has been run. This file states what was actually run; the win-builder line
above must be replaced with what those runs reported before the package is submitted, not with what
they would be expected to report.

## R CMD check results

Two notes locally.

* New submission. Expected and unavoidable for a first submission. The DESCRIPTION URLs previously
  returned 404 and no longer do; the repository now exists and both the homepage and the issues page
  resolve.
* "Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy" and
  "Skipping checking math rendering: package 'V8' unavailable". Both are missing local tools rather
  than properties of the package, and should not appear on CRAN's machines.

## Notes for the reviewer

The package compiles a small pure-C back-end with OpenMP. Each replicate draws from its own stream,
seeded from the base seed and the replicate index, so results do not depend on the number of
threads; there is a test asserting exactly that.

The correctness claim is exactness of the conditional reference distribution, so the tests do not
merely check that the code runs. The central test compares the sampler against the conditional law
obtained by COMPLETE ENUMERATION of all tables with the given margins, atom by atom, and a second
test compares a two by two case against `fisher.test`.
