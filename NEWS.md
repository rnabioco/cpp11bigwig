# cpp11bigwig (development version)

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
