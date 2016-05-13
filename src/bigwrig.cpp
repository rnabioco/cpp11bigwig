#include <Rcpp.h>
using namespace Rcpp ;

#include "bigWig.h"

//[[Rcpp::export]]
DataFrame read_bigwig_impl(std::string bwfname, std::string chrom, int start, int end) {

  //http://stackoverflow.com/questions/347949/how-to-convert-a-stdstring-to-const-char-or-char
  std::vector<char> bwfile(bwfname.begin(), bwfname.end()) ;
  bwfile.push_back('\0') ;

  const char mode = 'r' ;

  bigWigFile_t *bwf = NULL;
  // XXX change NULL to a CURL callback. see libBigWig demos
  bwf = bwOpen(&bwfile[0], NULL, &mode) ;

  if (!bwf)
    stop("Failed to open file: %s\n", bwfname) ;

  std::vector<std::string> chroms ;
  std::vector<int> starts ;
  std::vector<int> ends ;
  std::vector<float> vals;

  bwOverlappingIntervals_t *intervals = NULL ;

  int nchrom = bwf->cl->nKeys ;
  for (int nc = 0; nc<nchrom; ++nc) {

    char* cur_chrom = bwf->cl->chrom[nc] ;
    std::string cur_chrom_c = cur_chrom ;

    if (!chrom.empty() && chrom != cur_chrom_c) continue ;

    // set maximum boundaries of start / end are not specified
    if (start == -1) start = 0 ;
    if (end == -1) end = bwf->hdr->nBasesCovered ;

    intervals = bwGetValues(bwf, cur_chrom, start, end, 0) ;

    if (!intervals)
      stop("Failed to retreived intervals for %s\n", chrom) ;

    int nint = intervals->l ;
    for(int i=0; i<nint; ++i) {

      // +1 for 1-based coordinates
      int start = intervals->start[i] + 1;
      int end = start + 1 ;
      float val = intervals->value[i] ;

      starts.push_back(start) ;
      ends.push_back(end) ;
      vals.push_back(val) ;

      chroms.push_back(cur_chrom_c) ;
    }

    bwDestroyOverlappingIntervals(intervals) ;
  }

  return DataFrame::create( Named("chrom") = chroms,
                            Named("start") = starts,
                            Named("end") = ends,
                            Named("value") = vals) ;
}

/***R
library(dplyr)
read_bigwig_impl('src/libBigWig/test/test.bw', '', -1, -1) %>% as_data_frame
read_bigwig_impl('src/libBigWig/test/test.bw', '1', -1, -1) %>% as_data_frame
read_bigwig_impl('src/libBigWig/test/test.bw', 'x', -1, -1) %>% as_data_frame
*/
