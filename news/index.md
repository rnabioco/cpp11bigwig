# Changelog

## cpp11bigwig (development version)

- Fix a CRAN `gcc-ASAN` global-buffer-overflow reported when reading
  bigBed files. The autoSql schema parser no longer uses `std::regex`
  (which tripped an AddressSanitizer error inside libstdc++); it now
  parses the schema with simple string operations.

- [`read_bigwig()`](https://rnabioco.github.io/cpp11bigwig/reference/read_bigwig.md)
  and
  [`read_bigbed()`](https://rnabioco.github.io/cpp11bigwig/reference/read_bigbed.md)
  can now query multiple ranges in a single call. Pass equal-length (or
  length-1, recycled) `chrom`, `start`, and `end` vectors, or a
  `GRanges` of regions via `chrom`. For `read_bigwig(as = "Rle")`, a
  multi-range query returns a named `RleList` with one element per range
  ([\#18](https://github.com/rnabioco/cpp11bigwig/issues/18)).

## cpp11bigwig 0.2.0

CRAN release: 2026-06-18

- [`read_bigwig()`](https://rnabioco.github.io/cpp11bigwig/reference/read_bigwig.md)
  gains `as = "Rle"`, returning a per-base run-length-encoded vector
  spanning the queried range (an `Rle` for a single chromosome, or a
  named `RleList` for several). Uncovered bases are set to the `fill`
  value (default `0`; use `NA` to mark them missing)
  ([\#18](https://github.com/rnabioco/cpp11bigwig/issues/18)).

- Fix remote access to large bigWig/bigBed files. The HTTP `Range`
  header was not being set, so servers returned the entire file,
  crashing R or failing to open files larger than the read buffer
  ([\#18](https://github.com/rnabioco/cpp11bigwig/issues/18)).

## cpp11bigwig 0.1.3

CRAN release: 2025-12-11

- bigBed columns are now automatically coerced based on autoSql types
  ([\#12](https://github.com/rnabioco/cpp11bigwig/issues/12))

## cpp11bigwig 0.1.2

CRAN release: 2025-08-18

- Fix parsing of bigBed end coordinates.

## cpp11bigwig 0.1.1

CRAN release: 2025-01-19

- Sync with [libBigWig
  v0.4.8](https://github.com/dpryan79/libBigWig/releases/tag/0.4.8)

## cpp11bigwig 0.1.0

CRAN release: 2025-01-07

- Initial CRAN submission.

## cpp11bigwig 0.0.0.9000

- libBigWig is mostly untouched excpect for removal of `fprintf`
  statements (which R won’t allow in linked libraries) and fixups for
  ASAN errors, mostly GNU-specific pointer arithmetic. cpp11bigwig
  passes both ASAN and valgrind checks (via rhub).
