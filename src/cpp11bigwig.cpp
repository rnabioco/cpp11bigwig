#include <cpp11.hpp>
#include <cpp11/data_frame.hpp>

using namespace cpp11;

#include <string>
#include <vector>

#include "bigWig.h"

[[cpp11::register]]
writable::data_frame read_bigwig_cpp(std::string bwfname, sexp chrom, sexp start, sexp end, sexp raw) {

  //http://stackoverflow.com/questions/347949/how-to-convert-a-stdstring-to-const-char-or-char
  std::vector<char> bwfile(bwfname.begin(), bwfname.end()) ;
  bwfile.push_back('\0') ;

  const char mode = 'r' ;

  bigWigFile_t *bwf = NULL;

  // NULL can be a CURL callback. see libBigWig demos
  bwf = bwOpen(&bwfile[0], NULL, &mode) ;

  if (!bwf)
    stop("Failed to open file: '%s'\n", bwfname.c_str()) ;

  std::vector<std::string> chroms ;
  std::vector<int> starts ;
  std::vector<int> ends ;
  std::vector<float> vals;

  bwOverlappingIntervals_t *intervals = NULL ;

  int nchrom = bwf->cl->nKeys ;
  for (int nc = 0; nc<nchrom; ++nc) {

    char* bw_chrom = bwf->cl->chrom[nc] ;
    std::string bw_chrom_c(bw_chrom) ;

    if (!Rf_isNull(chrom)) {
      std::string r_chrom = as_cpp<std::string>(chrom) ;
      if (r_chrom != bw_chrom_c) continue ;
    }

    // set maximum boundaries if start / end are not specified
    int bw_start = Rf_isNull(start) ? 0 : as_cpp<int>(start) ;
    int bw_end = Rf_isNull(end) ? bwf->cl->len[nc] : as_cpp<int>(end) ;

    intervals = bwGetValues(bwf, bw_chrom, bw_start, bw_end, 0) ;

    if (!intervals)
      stop("Failed to retreived intervals for chrom `%s`\n", bw_chrom) ;

    int nint = intervals->l ;

    for(int i=0; i<nint; ++i) {

      int start = intervals->start[i] ;
      int end = start + 1 ;
      float val = intervals->value[i] ;

      starts.push_back(start) ;
      ends.push_back(end) ;
      vals.push_back(val) ;

      chroms.push_back(bw_chrom_c) ;
    }

    bwDestroyOverlappingIntervals(intervals) ;
  }

  bwClose(bwf) ;
  bwCleanup() ;

  if (raw) {
    return writable::data_frame({
        "chrom"_nm = chroms,
        "start"_nm = starts,
        "end"_nm = ends,
        "value"_nm = vals
      }) ;
  }

  std::vector<std::string> bg_chroms ;
  std::vector<int> bg_starts ;
  std::vector<int> bg_ends ;
  std::vector<float> bg_vals ;

  // collapse single-base intervals to interval spans with the same value
  int n = chroms.size() ;
  for (int i=0; i<n; ++i) {
    if (i == 0) {
      bg_chroms.push_back(chroms[i]) ;
      bg_starts.push_back(starts[i]) ;
      bg_ends.push_back(ends[i]) ;
      bg_vals.push_back(vals[i]) ;
    } else {
      if (chroms[i] == bg_chroms.back() && vals[i] == bg_vals.back()) {
        bg_ends.back() = ends[i] ;
      } else {
        bg_chroms.push_back(chroms[i]) ;
        bg_starts.push_back(starts[i]) ;
        bg_ends.push_back(ends[i]) ;
        bg_vals.push_back(vals[i]) ;
      }
    }
  }
  return writable::data_frame({
    "chrom"_nm = bg_chroms,
    "start"_nm = bg_starts,
    "end"_nm = bg_ends,
    "value"_nm = bg_vals
  }) ;
}
