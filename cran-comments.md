## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a minor release that fixes gcc-ASAN issues from the previous release,
  as well as a gcc-san (UBSan) "load of misaligned address" runtime error in the
  vendored libBigWig bigBed reader surfaced by the 0.3.0 pretest.

* The rchk issue in the cpp11bigwig CRAN checks will be resolved upstream in cpp11 <https://github.com/r-lib/cpp11/issues/408> and also affects other packages using cpp11 (e.g. cran.r-project.org/web/checks/check_results_tweenr.html).
