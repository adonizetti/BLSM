## Test environments
* local OS X install, R 3.5.0
* ubuntu 12.04 (on travis-ci), R 3.5.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.


## Resubmission

* Rephrased starting sentence in the description
* Added missing DOI in the description
* Added new examples and improved the existing ones. It's not possible for the main functions to give meaningful outputs if the simulations are too short, so I've introduced a couple of dontshow examples for testing purposes on a few iterations. I've also removed all the dontrun statements from examples which can be quickly executed. 

### Previous comments by CRAN reviewers

### Submission n. 2
Thanks, please do not start your description with "This package", package name or similar. Just start with "Computes..." or "Privides...".

Please add a DOI, arXiv or ISBN to your reference in your description text in the form
Hoff, Raftery and Handcock (2002) <doi:...>
Hoff, Raftery and Handcock (2002) <arXiv:...>
Hoff, Raftery and Handcock (2002, ISBN:...)
with no space after 'doi:', 'arXiv:' and angle brackets for auto-linking.

Most of your examples are wrapped in \dontrun{} and hence not tested. Please unwrap the examples if that is feasible and if they can be executed in < 5 sec for each Rd file or create additionally small toy examples. Something like
\examples{
       examples for users and checks:
       executable in < 5 sec
       \dontshow{
              examples for checks:
              executable in < 5 sec together with the examples above
              not shown to users
       }
       donttest{
              further examples for users (not used for checks)
       }
}
would be desirable.


### Submission n. 1 
The Description field contains

https://www.stat.washington.edu/people/pdhoff/Code/hoff_raftery_handcock_2002_jasa/.
Please enclose URLs in angle brackets (<...>).
