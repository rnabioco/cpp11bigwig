#include <cpp11.hpp>
#include <cpp11/data_frame.hpp>
using namespace cpp11;

#include "bigWig.h"

[[cpp11::register]]
writable::data_frame read_bigwig_impl(std::string bwfname, std::string chrom, int start, int end) {

  //http://stackoverflow.com/questions/347949/how-to-convert-a-stdstring-to-const-char-or-char
  std::vector<char> bwfile(bwfname.begin(), bwfname.end()) ;
  bwfile.push_back('\0') ;

  const char mode = 'r' ;

  bigWigFile_t *bwf = NULL;
  // XXX change NULL to a CURL callback. see libBigWig demos
  bwf = bwOpen(&bwfile[0], NULL, &mode) ;

  if (!bwf)
    stop("Failed to open file: %s\n", bwfname.c_str()) ;

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
      stop("Failed to retreived intervals for %s\n", chrom.c_str()) ;

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

  return writable::data_frame({
    "chrom"_nm = chroms,
    "start"_nm = starts,
    "end"_nm = ends,
    "value"_nm = vals
  }) ;

}
