# cpp11bigwig 0.3.0

* Fix a CRAN `gcc-san` (UBSan) `load of misaligned address` runtime error when
  reading a bigBed block that packs more than one record. In libBigWig's
  `bwValues.c`, records are stored as three `uint32_t` fields followed by a
  variable-length name, so every record after the first starts on an
  unaligned offset; the fields are now read with `memcpy` instead of an
  aligned `uint32_t` cast.

* Multi-range queries now open the file once per call instead of re-opening it
  for every range. The per-range loop moved into C++, so a query of many ranges
  (and especially a remote file, where each open re-fetches headers) is
  substantially faster.

* `read_bigbed()` no longer crashes on a bigBed file with no embedded autoSql
  schema. `bbGetSQL()` returns `NULL` in that case, and constructing a
  `std::string` from it was undefined behavior; such files now read back their
  `chrom`/`start`/`end` columns with no extra typed fields.

* The bigWig/bigBed readers now release the libBigWig file handle and read
  buffer when they error out (e.g. on an unreadable file or a failed interval
  query), rather than leaking them.

* Fix a CRAN `gcc-ASAN` global-buffer-overflow reported when reading bigBed
  files. The autoSql schema parser no longer uses `std::regex` (which tripped
  an AddressSanitizer error inside libstdc++); it now parses the schema with
  simple string operations.

* `read_bigwig()` and `read_bigbed()` can now query multiple ranges in a single
  call. Pass equal-length (or length-1, recycled) `chrom`, `start`, and `end`
  vectors, or a `GRanges` of regions via `chrom`. For `read_bigwig(as = "Rle")`,
  a multi-range query returns a named `RleList` with one element per range
  (#18).

# cpp11bigwig 0.2.0

* `read_bigwig()` gains `as = "Rle"`, returning a per-base run-length-encoded
  vector spanning the queried range (an `Rle` for a single chromosome, or a
  named `RleList` for several). Uncovered bases are set to the `fill` value
  (default `0`; use `NA` to mark them missing) (#18).

* Fix remote access to large bigWig/bigBed files. The HTTP `Range` header was
  not being set, so servers returned the entire file, crashing R or failing to
  open files larger than the read buffer (#18).

# cpp11bigwig 0.1.3

* bigBed columns are now automatically coerced based on autoSql types (#12)

# cpp11bigwig 0.1.2

* Fix parsing of bigBed end coordinates.
  
# cpp11bigwig 0.1.1

* Sync with [libBigWig v0.4.8](https://github.com/dpryan79/libBigWig/releases/tag/0.4.8)

# cpp11bigwig 0.1.0

* Initial CRAN submission.

# cpp11bigwig 0.0.0.9000

* libBigWig is mostly untouched excpect for removal of `fprintf` statements (which R won't allow in linked libraries) and fixups for ASAN errors, mostly GNU-specific pointer arithmetic. cpp11bigwig passes both ASAN and valgrind checks (via rhub).
