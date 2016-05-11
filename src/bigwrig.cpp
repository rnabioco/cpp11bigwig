#include <Rcpp.h>
using namespace Rcpp ;

#include "bigWig.h"

// [[Rcpp::export]]
DataFrame read_bigwig_impl(std::string bwfname, std::string chrom, int start, int end) {

  //http://stackoverflow.com/questions/347949/how-to-convert-a-stdstring-to-const-char-or-char
  std::vector<char> bwfile(bwfname.begin(), bwfname.end()) ;
  bwfile.push_back('\0') ;

  const char mode = 'r' ;

  bigWigFile_t *bwf = NULL;
  bwf = bwOpen(&bwfile[0], NULL, &mode) ;

  if (!bwf)
    stop("Failed to open file: %s\n", bwfname) ;

  std::vector<std::string> chroms ;
  std::vector<int> starts ;
  std::vector<int> ends ;

  // chrom list
  int nkeys = bwf->cl->nKeys ;

  for(int i=0; i<nkeys; ++i) {

    bwOverlappingIntervals_t *intervals = NULL ;
    intervals =  bwGetValues(bwf, bwf->cl->chrom[i], start, end, 0) ;

    if (!intervals)
      stop("Failed to retreived intervals for %s", bwf->cl->chrom[i]) ;

    int nint = intervals->l ;

    for(int i; i<nint; ++i) {
      int start = intervals->start[i] ;
      starts[i] = start ;
      ends[i] = start + 1 ;
    }

    bwDestroyOverlappingIntervals(intervals) ;
  }

  return DataFrame::create( Named("chrom") = chroms,
                            Named("start") = starts,
                            Named("end") = ends) ;
}

/***R

*/
