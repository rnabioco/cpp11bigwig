#include <cpp11.hpp>
#include <cpp11/data_frame.hpp>

using namespace cpp11;

#include <cmath>
#include <sstream>
#include <vector>

#include "libBigWig/bigWig.h"
// bwCommon.h (bwSetPos/bwRead, used to verify file magic) lacks extern "C"
// guards, so declare its symbols with C linkage to match the C library.
extern "C" {
#include "libBigWig/bwCommon.h"
}

// R type categories for autoSql types
enum class RType { Integer, Double, String };

// Determine R type from autoSql type string
RType autosql_to_rtype(const std::string& type) {
  // Strip array notation like int[blockCount]
  std::string base_type = type;
  size_t bracket_pos = type.find('[');
  if (bracket_pos != std::string::npos) {
    // Array types stay as strings (comma-separated values)
    return RType::String;
  }

  if (type == "uint" || type == "int") {
    return RType::Integer;
  } else if (type == "float" || type == "double") {
    return RType::Double;
  }
  return RType::String;
}

// Parse autoSql schema to extract field names and types
// Returns vector of pairs: (name, type)
std::vector<std::pair<std::string, std::string>> parse_autosql(const std::string& sql) {
  std::vector<std::pair<std::string, std::string>> fields;

  std::istringstream stream(sql);
  std::string line;

  // Match lines like: "   uint   chromStart;  ..." or "   string name;  ..."
  // (formerly a std::regex `^\s*(\S+)\s+(\S+);` — parsed by hand here to avoid
  // a libstdc++ std::regex global-buffer-overflow flagged by AddressSanitizer).
  while (std::getline(stream, line)) {
    std::istringstream ls(line);
    std::string type, name_tok;

    // need at least two whitespace-delimited tokens: the type and the name
    if (!(ls >> type >> name_tok))
      continue;

    // the name token must end in a ';' (matches the trailing `;` in the regex);
    // use the last ';' to mirror the greedy `\S+;`
    size_t semi = name_tok.rfind(';');
    if (semi == std::string::npos)
      continue;

    std::string name = name_tok.substr(0, semi);
    if (name.empty())
      continue;

    fields.push_back({name, type});
  }

  return fields;
}

// Canonical UCSC BED fields (chrom, chromStart, chromEnd, then the bed4..bed12
// columns), used as a fallback when a bigBed has no embedded autoSql schema.
// Mirrors the standard bed12 autoSql so a schema-less file reads the same as one
// encoded with `-as=bed12.as`. This naming is a cpp11bigwig convenience: UCSC's
// `bigBedToBed -header`/`-tsv` instead aborts when a file has no autoSql rather
// than inferring names, so read_bigbed() emits a message in that case. Sized to
// the file header's field counts:
// `field_count` is the total number of columns and `defined_field_count` is how
// many of the fixed-format BED fields are present (3-12); any field beyond that
// is a bedN+ extra column of unknown type and falls back to a generic string.
std::vector<std::pair<std::string, std::string>> default_bed_fields(uint16_t field_count,
                                                                    uint16_t defined_field_count) {
  static const std::pair<const char*, const char*> bed_fields[] = {
      {"chrom", "string"},
      {"start", "uint"},
      {"end", "uint"},
      {"name", "string"},
      {"score", "uint"},
      {"strand", "char[1]"},
      {"thickStart", "uint"},
      {"thickEnd", "uint"},
      {"itemRgb", "uint"},
      {"blockCount", "int"},
      {"blockSizes", "int[blockCount]"},
      {"blockStarts", "int[blockCount]"}};

  std::vector<std::pair<std::string, std::string>> fields;
  for (uint16_t i = 0; i < field_count; ++i) {
    if (i < defined_field_count && i < 12) {
      fields.push_back({bed_fields[i].first, bed_fields[i].second});
    } else {
      // bedN+ extra field with unknown type: keep as string, generic name
      fields.push_back({"field" + std::to_string(i + 1), "string"});
    }
  }
  return fields;
}

// Split a string by delimiter
std::vector<std::string> split_string(const std::string& str, char delim) {
  std::vector<std::string> tokens;
  std::istringstream stream(str);
  std::string token;
  while (std::getline(stream, token, delim)) {
    tokens.push_back(token);
  }
  return tokens;
}

// Each query is one (chrom, start, end) range; the inputs are equal-length
// vectors describing several ranges. An NA `chrom` selects every chromosome;
// an NA `start`/`end` defaults to the chromosome's bounds. The file is opened
// once and every range is served from that handle. Returns a list with one
// data frame per range (the R layer combines or post-processes them).
[[cpp11::register]]
writable::list read_bigwig_cpp(std::string bwfname, strings chroms_r, integers starts_r,
                               integers ends_r) {
  const char* bwfile = bwfname.c_str();

  bigWigFile_t* bwf = NULL;

  // initialize libBigWig (allocates the read buffer required for remote files)
  if (bwInit(1 << 17) != 0)
    stop("Failed to initialize libBigWig\n");

  // 2nd arg is an optional curl-options callback; NULL uses libBigWig's
  // defaults. Remote (http/https/ftp) URLs still work — see libBigWig demos.
  bwf = bwOpen(bwfile, NULL, "r");

  if (!bwf) {
    bwCleanup();
    stop("Failed to open file: '%s'\n", bwfname.c_str());
  }

  R_xlen_t nranges = chroms_r.size();
  writable::list out(nranges);

  int nchrom = bwf->cl->nKeys;
  for (R_xlen_t r = 0; r < nranges; ++r) {
    bool any_chrom = is_na(chroms_r[r]);
    std::string want_chrom = any_chrom ? std::string() : static_cast<std::string>(chroms_r[r]);
    bool start_na = starts_r[r] == NA_INTEGER;
    bool end_na = ends_r[r] == NA_INTEGER;

    std::vector<std::string> chroms;
    std::vector<int> starts;
    std::vector<int> ends;
    std::vector<float> vals;

    for (int nc = 0; nc < nchrom; ++nc) {
      char* bw_chrom = bwf->cl->chrom[nc];
      std::string bw_chrom_c(bw_chrom);

      if (!any_chrom && want_chrom != bw_chrom_c)
        continue;

      // set maximum boundaries if start / end are not specified
      int bw_start = start_na ? 0 : starts_r[r];
      int bw_end = end_na ? bwf->cl->len[nc] : ends_r[r];

      bwOverlappingIntervals_t* intervals = bwGetValues(bwf, bw_chrom, bw_start, bw_end, 0);

      if (!intervals) {
        bwClose(bwf);
        bwCleanup();
        stop("Failed to retrieve intervals for chrom `%s`\n", bw_chrom);
      }

      int nint = intervals->l;

      for (int i = 0; i < nint; ++i) {
        int start = intervals->start[i];
        int end = start + 1;
        float val = intervals->value[i];

        if (i == 0) {
          chroms.push_back(bw_chrom_c);
          starts.push_back(start);
          ends.push_back(end);
          vals.push_back(val);
        } else {
          if (start == ends.back() && val == vals.back()) {
            ends.back() = end;
          } else {
            chroms.push_back(bw_chrom_c);
            starts.push_back(start);
            ends.push_back(end);
            vals.push_back(val);
          }
        }
      }

      bwDestroyOverlappingIntervals(intervals);
    }

    out[r] = writable::data_frame(
        {"chrom"_nm = chroms, "start"_nm = starts, "end"_nm = ends, "value"_nm = vals});
  }

  bwClose(bwf);
  bwCleanup();

  return out;
}

// As with read_bigwig_cpp, the inputs are equal-length per-range vectors and
// the file is opened once. The autoSql schema (and thus the column layout) is
// parsed a single time up front and shared across ranges. Returns a list with
// one data frame per range.
[[cpp11::register]]
writable::list read_bigbed_cpp(std::string bbfname, strings chroms_r, integers starts_r,
                               integers ends_r) {
  const char* bbfile = bbfname.c_str();
  bigWigFile_t* bbf = NULL;

  // initialize libBigWig (allocates the read buffer required for remote files)
  if (bwInit(1 << 17) != 0)
    stop("Failed to initialize libBigWig\n");

  // 2nd arg is an optional curl-options callback; NULL uses libBigWig's
  // defaults. Remote (http/https/ftp) URLs still work — see libBigWig demos.
  bbf = bbOpen(bbfile, NULL);

  if (!bbf) {
    bwCleanup();
    stop("Failed to open file: '%s'\n", bbfname.c_str());
  }

  // Get autoSql schema and parse field info. bbGetSQL() returns NULL for
  // bigBed files without an embedded schema (or on read error); guard against
  // constructing a std::string from NULL.
  char* sql = bbGetSQL(bbf);
  std::string sql_str(sql ? sql : "");
  free(sql);

  auto fields = parse_autosql(sql_str);

  // Without an embedded schema (or on a partial parse) the field list comes up
  // short of the layout recorded in the file header. Fall back to the canonical
  // BED column names/types sized to the header counts so a schema-less bed12
  // still yields all 12 columns instead of just chrom/start/end. A complete
  // schema parses to exactly `fieldCount` entries and is never overridden.
  if (fields.size() < bbf->hdr->fieldCount) {
    fields = default_bed_fields(bbf->hdr->fieldCount, bbf->hdr->definedFieldCount);
  }

  // Skip first 3 fields (chrom, chromStart, chromEnd) - handled separately
  size_t num_extra_fields = fields.size() > 3 ? fields.size() - 3 : 0;

  // Determine R types for extra fields
  std::vector<RType> field_rtypes;
  std::vector<std::string> field_names;
  for (size_t i = 3; i < fields.size(); ++i) {
    field_names.push_back(fields[i].first);
    field_rtypes.push_back(autosql_to_rtype(fields[i].second));
  }

  R_xlen_t nranges = chroms_r.size();
  writable::list out(nranges);

  int nchrom = bbf->cl->nKeys;
  for (R_xlen_t r = 0; r < nranges; ++r) {
    bool any_chrom = is_na(chroms_r[r]);
    std::string want_chrom = any_chrom ? std::string() : static_cast<std::string>(chroms_r[r]);
    bool start_na = starts_r[r] == NA_INTEGER;
    bool end_na = ends_r[r] == NA_INTEGER;

    // Storage for coordinate columns
    std::vector<std::string> chroms;
    std::vector<int> starts;
    std::vector<int> ends;

    // Storage for typed extra columns
    std::vector<std::vector<int>> int_cols(num_extra_fields);
    std::vector<std::vector<double>> dbl_cols(num_extra_fields);
    std::vector<std::vector<std::string>> str_cols(num_extra_fields);

    for (int nc = 0; nc < nchrom; ++nc) {
      char* bb_chrom = bbf->cl->chrom[nc];
      std::string bb_chrom_c(bb_chrom);

      if (!any_chrom && want_chrom != bb_chrom_c)
        continue;

      // set maximum boundaries if start / end are not specified
      int bb_start = start_na ? 0 : starts_r[r];
      int bb_end = end_na ? bbf->cl->len[nc] : ends_r[r];

      bbOverlappingEntries_t* intervals =
          bbGetOverlappingEntries(bbf, bb_chrom, bb_start, bb_end, 1);

      if (!intervals) {
        bwClose(bbf);
        bwCleanup();
        stop("Failed to retrieve intervals for chrom `%s`\n", bb_chrom);
      }

      int nint = intervals->l;

      for (int i = 0; i < nint; ++i) {
        int ivl_start = intervals->start[i];
        int ivl_end = intervals->end[i];
        std::string val_str = intervals->str[i];

        chroms.push_back(bb_chrom_c);
        starts.push_back(ivl_start);
        ends.push_back(ivl_end);

        // Parse tab-separated extra fields
        std::vector<std::string> tokens = split_string(val_str, '\t');

        for (size_t j = 0; j < num_extra_fields; ++j) {
          std::string token = (j < tokens.size()) ? tokens[j] : "";

          switch (field_rtypes[j]) {
            case RType::Integer:
              if (token.empty()) {
                int_cols[j].push_back(NA_INTEGER);
              } else {
                try {
                  int_cols[j].push_back(std::stoi(token));
                } catch (...) {
                  int_cols[j].push_back(NA_INTEGER);
                }
              }
              break;
            case RType::Double:
              if (token.empty()) {
                dbl_cols[j].push_back(NA_REAL);
              } else {
                try {
                  dbl_cols[j].push_back(std::stod(token));
                } catch (...) {
                  dbl_cols[j].push_back(NA_REAL);
                }
              }
              break;
            case RType::String:
              str_cols[j].push_back(token);
              break;
          }
        }
      }

      bbDestroyOverlappingEntries(intervals);
    }

    // Build the data frame with named columns
    writable::list result;
    writable::strings col_names;

    // Add coordinate columns
    result.push_back(writable::strings(chroms.begin(), chroms.end()));
    col_names.push_back("chrom");

    result.push_back(writable::integers(starts.begin(), starts.end()));
    col_names.push_back("start");

    result.push_back(writable::integers(ends.begin(), ends.end()));
    col_names.push_back("end");

    // Add typed extra columns
    for (size_t j = 0; j < num_extra_fields; ++j) {
      switch (field_rtypes[j]) {
        case RType::Integer:
          result.push_back(writable::integers(int_cols[j].begin(), int_cols[j].end()));
          break;
        case RType::Double:
          result.push_back(writable::doubles(dbl_cols[j].begin(), dbl_cols[j].end()));
          break;
        case RType::String:
          result.push_back(writable::strings(str_cols[j].begin(), str_cols[j].end()));
          break;
      }
      col_names.push_back(field_names[j]);
    }

    result.attr("names") = col_names;
    result.attr("class") = writable::strings({"data.frame"});
    result.attr("row.names") = writable::integers({NA_INTEGER, -static_cast<int>(chroms.size())});

    out[r] = as_cpp<writable::data_frame>(result);
  }

  bwClose(bbf);
  bwCleanup();

  // Signal whether the file carried an embedded autoSql schema. sql_str is ""
  // exactly when bbGetSQL() returned NULL (no schema), in which case the column
  // names above came from the default_bed_fields() fallback rather than the
  // file. read_bigbed() reads this off the returned list to emit a one-time
  // message; it rides as an attribute so the C++ signature stays unchanged.
  out.attr("has_autosql") = writable::logicals({r_bool(!sql_str.empty())});

  return out;
}

[[cpp11::register]]
std::string bigbed_sql_cpp(std::string bbfname) {
  const char* bbfile = bbfname.c_str();

  bigWigFile_t* bbf = NULL;

  // initialize libBigWig (allocates the read buffer required for remote files)
  if (bwInit(1 << 17) != 0)
    stop("Failed to initialize libBigWig\n");

  // 2nd arg is an optional curl-options callback; NULL uses libBigWig's
  // defaults. Remote (http/https/ftp) URLs still work — see libBigWig demos.
  bbf = bbOpen(bbfile, NULL);

  if (!bbf) {
    bwCleanup();
    stop("Failed to open file: '%s'\n", bbfname.c_str());
  }

  char* sql = bbGetSQL(bbf);

  // NULL when the bigBed has no embedded autoSql schema (or on read error)
  std::string sql_str(sql ? sql : "");
  free(sql);

  bwClose(bbf);
  bwCleanup();

  return (sql_str);
}

// Report header metadata for a bigBed file. `field_count` is the total number of
// columns and `defined_field_count` is how many of the fixed-format BED columns
// are present (3-12) -- the authoritative way to identify the BED variant (e.g.
// `defined_field_count == 12` is a genuine BED12). `autosql` is the embedded
// schema string (empty when the file has none).
[[cpp11::register]]
writable::list bigbed_info_cpp(std::string bbfname) {
  const char* bbfile = bbfname.c_str();

  // initialize libBigWig (allocates the read buffer required for remote files)
  if (bwInit(1 << 17) != 0)
    stop("Failed to initialize libBigWig\n");

  // 2nd arg is an optional curl-options callback; NULL uses libBigWig's
  // defaults. Remote (http/https/ftp) URLs still work — see libBigWig demos.
  bigWigFile_t* bbf = bbOpen(bbfile, NULL);

  if (!bbf) {
    bwCleanup();
    stop("Failed to open file: '%s'\n", bbfname.c_str());
  }

  // bbOpen() accepts either magic, so verify the file is really a bigBed by
  // re-reading the magic from the open handle (avoids a second remote open)
  uint32_t magic = 0;
  if (bwSetPos(bbf, 0) != 0 || bwRead(&magic, sizeof(uint32_t), 1, bbf) != 1 ||
      magic != BIGBED_MAGIC) {
    bwClose(bbf);
    bwCleanup();
    stop("Not a bigBed file: '%s'\n", bbfname.c_str());
  }

  char* sql = bbGetSQL(bbf);
  std::string sql_str(sql ? sql : "");
  free(sql);

  writable::list out;
  writable::strings nms;

  out.push_back(writable::integers({static_cast<int>(bbf->hdr->version)}));
  nms.push_back("version");
  out.push_back(writable::integers({static_cast<int>(bbf->cl->nKeys)}));
  nms.push_back("n_chroms");
  out.push_back(writable::integers({static_cast<int>(bbf->hdr->fieldCount)}));
  nms.push_back("field_count");
  out.push_back(writable::integers({static_cast<int>(bbf->hdr->definedFieldCount)}));
  nms.push_back("defined_field_count");
  out.push_back(writable::doubles({static_cast<double>(bbf->hdr->nBasesCovered)}));
  nms.push_back("n_bases_covered");
  out.push_back(writable::strings({sql_str}));
  nms.push_back("autosql");

  out.attr("names") = nms;

  bwClose(bbf);
  bwCleanup();

  return out;
}

// Report header metadata and summary statistics for a bigWig file.
[[cpp11::register]]
writable::list bigwig_info_cpp(std::string bwfname) {
  const char* bwfile = bwfname.c_str();

  // initialize libBigWig (allocates the read buffer required for remote files)
  if (bwInit(1 << 17) != 0)
    stop("Failed to initialize libBigWig\n");

  // 2nd arg is an optional curl-options callback; NULL uses libBigWig's
  // defaults. Remote (http/https/ftp) URLs still work — see libBigWig demos.
  bigWigFile_t* bw = bwOpen(bwfile, NULL, "r");

  if (!bw) {
    bwCleanup();
    stop("Failed to open file: '%s'\n", bwfname.c_str());
  }

  // bwOpen() accepts either magic, so verify the file is really a bigWig by
  // re-reading the magic from the open handle (avoids a second remote open)
  uint32_t magic = 0;
  if (bwSetPos(bw, 0) != 0 || bwRead(&magic, sizeof(uint32_t), 1, bw) != 1 ||
      magic != BIGWIG_MAGIC) {
    bwClose(bw);
    bwCleanup();
    stop("Not a bigWig file: '%s'\n", bwfname.c_str());
  }

  // file-level summary stats live in the header; derive mean/std from the sums
  double nbases = static_cast<double>(bw->hdr->nBasesCovered);
  double mean = nbases > 0 ? bw->hdr->sumData / nbases : NA_REAL;
  double sd = NA_REAL;
  if (nbases > 0) {
    double var = (bw->hdr->sumSquared - bw->hdr->sumData * bw->hdr->sumData / nbases) / nbases;
    sd = var > 0 ? std::sqrt(var) : 0.0;
  }

  writable::list out;
  writable::strings nms;

  out.push_back(writable::integers({static_cast<int>(bw->hdr->version)}));
  nms.push_back("version");
  out.push_back(writable::integers({static_cast<int>(bw->hdr->nLevels)}));
  nms.push_back("n_levels");
  out.push_back(writable::integers({static_cast<int>(bw->cl->nKeys)}));
  nms.push_back("n_chroms");
  out.push_back(writable::doubles({nbases}));
  nms.push_back("n_bases_covered");
  out.push_back(writable::doubles({bw->hdr->minVal}));
  nms.push_back("min");
  out.push_back(writable::doubles({bw->hdr->maxVal}));
  nms.push_back("max");
  out.push_back(writable::doubles({mean}));
  nms.push_back("mean");
  out.push_back(writable::doubles({sd}));
  nms.push_back("std");

  out.attr("names") = nms;

  bwClose(bw);
  bwCleanup();

  return out;
}

// Report whether this build was compiled with libcurl (remote file) support.
[[cpp11::register]]
bool bigwig_has_curl_cpp() {
  return LIBBIGWIG_CURL == 1;
}
