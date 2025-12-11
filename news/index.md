# Changelog

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
