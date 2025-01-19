## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a patch release to sync with libBigWig v.0.4.8 <https://github.com/dpryan79/libBigWig/releases/tag/0.4.8>

* The rchk issue in the cpp11bigwig CRAN checks will be resolved upstream in cpp11 <https://github.com/r-lib/cpp11/issues/408> and also affects other packages using cpp11 (e.g. cran.r-project.org/web/checks/check_results_tweenr.html).

## ASAN / valgrind

cpp11bigwig passes [ASAN and valgrind checks](https://github.com/rnabioco/cpp11bigwig/actions/runs/12857168623).
