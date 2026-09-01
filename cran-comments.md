## Test environments

* local macOS 15 (arm64), R 4.5.1, `R CMD check --as-cran`: 0 errors, 0 warnings, 2 notes.
* win-builder R-devel (R Under development, 2026-08-31 r90457 ucrt): 0 errors, 0 warnings, 1 note.
  The note was "New submission" together with a spell-check flag on the word "conditionings", which
  has since been reworded. Notably the local "HTML version of manual" note did NOT appear there,
  confirming it was a missing local tool rather than a property of the package.
* win-builder R-release: submitted at the same time.

No GitHub Actions matrix has been run.

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
