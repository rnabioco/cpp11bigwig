#include <cpp11.hpp>
#include <cpp11/data_frame.hpp>

using namespace cpp11;

#include <vector>

#include "bigWig.h"

[[cpp11::register]]
writable::data_frame read_bigwig_cpp(std::string bwfname, sexp chrom, sexp start, sexp end) {

  const char* bwfile = bwfname.c_str() ;

  bigWigFile_t *bwf = NULL;

  // NULL can be a CURL callback. see libBigWig demos
  bwf = bwOpen(bwfile, NULL, "r") ;

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

      if (i == 0) {
        chroms.push_back(bw_chrom_c) ;
        starts.push_back(start) ;
        ends.push_back(end) ;
        vals.push_back(val) ;
      } else {
        if (start == ends.back() && val == vals.back()) {
          ends.back() = end ;
        } else {
          chroms.push_back(bw_chrom_c) ;
          starts.push_back(start) ;
          ends.push_back(end) ;
          vals.push_back(val) ;
        }
      }
    }

    bwDestroyOverlappingIntervals(intervals) ;
  }

  bwClose(bwf) ;
  bwCleanup() ;

  return writable::data_frame({
    "chrom"_nm = chroms,
    "start"_nm = starts,
    "end"_nm = ends,
    "value"_nm = vals
  }) ;
}

[[cpp11::register]]
writable::data_frame read_bigbed_cpp(std::string bbfname, sexp chrom, sexp start, sexp end) {

  const char* bbfile = bbfname.c_str() ;
  bigWigFile_t *bbf = NULL;

  // NULL can be a CURL callback. see libBigWig demos
  bbf = bbOpen(bbfile, NULL) ;

  if (!bbf)
    stop("Failed to open file: '%s'\n", bbfname.c_str()) ;

  std::vector<std::string> chroms ;
  std::vector<int> starts ;
  std::vector<int> ends ;
  std::vector<std::string> vals;

  bbOverlappingEntries_t *intervals ;

  int nchrom = bbf->cl->nKeys ;
  for (int nc = 0; nc<nchrom; ++nc) {

    char* bb_chrom = bbf->cl->chrom[nc] ;
    std::string bb_chrom_c(bb_chrom) ;

    if (!Rf_isNull(chrom)) {
      std::string r_chrom = as_cpp<std::string>(chrom) ;
      if (r_chrom != bb_chrom_c) continue ;
    }

    // set maximum boundaries if start / end are not specified
    int bb_start = Rf_isNull(start) ? 0 : as_cpp<int>(start) ;
    int bb_end = Rf_isNull(end) ? bbf->cl->len[nc] : as_cpp<int>(end) ;

    intervals = bbGetOverlappingEntries(bbf, bb_chrom, bb_start, bb_end, 1) ;

    if (!intervals)
      stop("Failed to retreived intervals for chrom `%s`\n", bb_chrom) ;

    int nint = intervals->l ;

    for(int i=0; i<nint; ++i) {

      int start = intervals->start[i] ;
      int end = start + 1 ;
      std::string val = intervals->str[i] ;

      if (i == 0) {
        chroms.push_back(bb_chrom_c) ;
        starts.push_back(start) ;
        ends.push_back(end) ;
        vals.push_back(val) ;
      } else {
        if (start == ends.back() && val == vals.back()) {
          ends.back() = end ;
        } else {
          chroms.push_back(bb_chrom_c) ;
          starts.push_back(start) ;
          ends.push_back(end) ;
          vals.push_back(val) ;
        }
      }
    }

    bbDestroyOverlappingEntries(intervals) ;
  }

  bwClose(bbf) ;
  bwCleanup() ;

  return writable::data_frame({
    "chrom"_nm = chroms,
      "start"_nm = starts,
      "end"_nm = ends,
      "value"_nm = vals
  }) ;
}

[[cpp11::register]]
std::string bigbed_sql_cpp(std::string bbfname) {

  const char* bbfile = bbfname.c_str() ;

  bigWigFile_t *bbf = NULL;

  // NULL can be a CURL callback. see libBigWig demos
  bbf = bwOpen(bbfile, NULL, "r") ;

  if (!bbf)
    stop("Failed to open file: '%s'\n", bbfname.c_str()) ;

  char* sql = bbGetSQL(bbf) ;

  std::string sql_str(sql) ;
  free(sql) ;

  bwClose(bbf) ;

  return(sql_str) ;
}
