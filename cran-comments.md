## Test environments
* local OS X install, R 3.4.4
* ubuntu 12.04 (on travis-ci), R 3.4.4
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

Resubmitting due to:
* missing DOI in the description - added
* too many dontrun examples. It's not possible for the main functions to give meaningful outputs if the simulations are too short, so I've added at least one dontshow example for testing purposes on a few iterations. I've also removed all the "\dontrun" from examples which should be quickly executable.