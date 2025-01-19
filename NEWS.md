# cpp11bigwig 0.1.1

* Sync with [libBigWig v0.4.8](https://github.com/dpryan79/libBigWig/releases/tag/0.4.8)

# cpp11bigwig 0.1.0

* Initial CRAN submission.

# cpp11bigwig 0.0.0.9000

* libBigWig is mostly untouched excpect for removal of `fprintf` statements (which R won't allow in linked libraries) and fixups for ASAN errors, mostly GNU-specific pointer arithmetic. cpp11bigwig passes both ASAN and valgrind checks (via rhub).
