#include <cpp11.hpp>
#include <cpp11/data_frame.hpp>

using namespace cpp11;

#include <regex>
#include <sstream>
#include <vector>

#include "libBigWig/bigWig.h"

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
  std::regex field_regex(R"(^\s*(\S+)\s+(\S+);)");

  while (std::getline(stream, line)) {
    std::smatch match;
    if (std::regex_search(line, match, field_regex)) {
      std::string type = match[1].str();
      std::string name = match[2].str();
      fields.push_back({name, type});
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

[[cpp11::register]]
writable::data_frame read_bigwig_cpp(std::string bwfname, sexp chrom, sexp start, sexp end) {
  const char* bwfile = bwfname.c_str();

  bigWigFile_t* bwf = NULL;

  // NULL can be a CURL callback. see libBigWig demos
  bwf = bwOpen(bwfile, NULL, "r");

  if (!bwf)
    stop("Failed to open file: '%s'\n", bwfname.c_str());

  std::vector<std::string> chroms;
  std::vector<int> starts;
  std::vector<int> ends;
  std::vector<float> vals;

  bwOverlappingIntervals_t* intervals = NULL;

  int nchrom = bwf->cl->nKeys;
  for (int nc = 0; nc < nchrom; ++nc) {
    char* bw_chrom = bwf->cl->chrom[nc];
    std::string bw_chrom_c(bw_chrom);

    if (!Rf_isNull(chrom)) {
      std::string r_chrom = as_cpp<std::string>(chrom);
      if (r_chrom != bw_chrom_c)
        continue;
    }

    // set maximum boundaries if start / end are not specified
    int bw_start = Rf_isNull(start) ? 0 : as_cpp<int>(start);
    int bw_end = Rf_isNull(end) ? bwf->cl->len[nc] : as_cpp<int>(end);

    intervals = bwGetValues(bwf, bw_chrom, bw_start, bw_end, 0);

    if (!intervals)
      stop("Failed to retreived intervals for chrom `%s`\n", bw_chrom);

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

  bwClose(bwf);
  bwCleanup();

  return writable::data_frame(
      {"chrom"_nm = chroms, "start"_nm = starts, "end"_nm = ends, "value"_nm = vals});
}

[[cpp11::register]]
writable::data_frame read_bigbed_cpp(std::string bbfname, sexp chrom, sexp start, sexp end) {
  const char* bbfile = bbfname.c_str();
  bigWigFile_t* bbf = NULL;

  // NULL can be a CURL callback. see libBigWig demos
  bbf = bbOpen(bbfile, NULL);

  if (!bbf)
    stop("Failed to open file: '%s'\n", bbfname.c_str());

  // Get autoSql schema and parse field info
  char* sql = bbGetSQL(bbf);
  std::string sql_str(sql);
  free(sql);

  auto fields = parse_autosql(sql_str);

  // Skip first 3 fields (chrom, chromStart, chromEnd) - handled separately
  size_t num_extra_fields = fields.size() > 3 ? fields.size() - 3 : 0;

  // Determine R types for extra fields
  std::vector<RType> field_rtypes;
  std::vector<std::string> field_names;
  for (size_t i = 3; i < fields.size(); ++i) {
    field_names.push_back(fields[i].first);
    field_rtypes.push_back(autosql_to_rtype(fields[i].second));
  }

  // Storage for coordinate columns
  std::vector<std::string> chroms;
  std::vector<int> starts;
  std::vector<int> ends;

  // Storage for typed extra columns
  std::vector<std::vector<int>> int_cols(num_extra_fields);
  std::vector<std::vector<double>> dbl_cols(num_extra_fields);
  std::vector<std::vector<std::string>> str_cols(num_extra_fields);

  bbOverlappingEntries_t* intervals;

  int nchrom = bbf->cl->nKeys;
  for (int nc = 0; nc < nchrom; ++nc) {
    char* bb_chrom = bbf->cl->chrom[nc];
    std::string bb_chrom_c(bb_chrom);

    if (!Rf_isNull(chrom)) {
      std::string r_chrom = as_cpp<std::string>(chrom);
      if (r_chrom != bb_chrom_c)
        continue;
    }

    // set maximum boundaries if start / end are not specified
    int bb_start = Rf_isNull(start) ? 0 : as_cpp<int>(start);
    int bb_end = Rf_isNull(end) ? bbf->cl->len[nc] : as_cpp<int>(end);

    intervals = bbGetOverlappingEntries(bbf, bb_chrom, bb_start, bb_end, 1);

    if (!intervals)
      stop("Failed to retrieve intervals for chrom `%s`\n", bb_chrom);

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

  bwClose(bbf);
  bwCleanup();

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

  return as_cpp<writable::data_frame>(result);
}

[[cpp11::register]]
std::string bigbed_sql_cpp(std::string bbfname) {
  const char* bbfile = bbfname.c_str();

  bigWigFile_t* bbf = NULL;

  // NULL can be a CURL callback. see libBigWig demos
  bbf = bwOpen(bbfile, NULL, "r");

  if (!bbf)
    stop("Failed to open file: '%s'\n", bbfname.c_str());

  char* sql = bbGetSQL(bbf);

  std::string sql_str(sql);
  free(sql);

  bwClose(bbf);

  return (sql_str);
}
