## Golay correction

## Nucleotide->bit mappings (decode allows N; encode does not)
DECODE_GOLAY_NT_TO_BITS <- c(A="11", C="00", T="10", G="01", N="11")
ENCODE_GOLAY_NT_TO_BITS <- c(A="11", C="00", T="10", G="01")

## Parity submatrix P (12x12), as in the original code
DEFAULT_P <- matrix(
  c(
    0,1,1,1,1,1,1,1,1,1,1,1,
    1,1,1,0,1,1,1,0,0,0,1,0,
    1,1,0,1,1,1,0,0,0,1,0,1,
    1,0,1,1,1,0,0,0,1,0,1,1,
    1,1,1,1,0,0,0,1,0,1,1,0,
    1,1,1,0,0,0,1,0,1,1,0,1,
    1,1,0,0,0,1,0,1,1,0,1,1,
    1,0,0,0,1,0,1,1,0,1,1,1,
    1,0,0,1,0,1,1,0,1,1,1,0,
    1,0,1,0,1,1,0,1,1,1,0,0,
    1,1,0,1,1,0,1,1,1,0,0,0,
    1,0,1,1,0,1,1,1,0,0,0,1
  ),
  nrow = 12, byrow = TRUE
)

## Generator (G) and parity-check (H) matrices
DEFAULT_G <- cbind(DEFAULT_P, diag(12))
DEFAULT_H <- cbind(diag(12), t(DEFAULT_P))

## Build all error vectors with Hamming weight <= 3 (length 24)
.make_3bit_errors <- function(veclen = 24L) {
  errs <- list()
  idx <- 1L
  z <- rep(0L, veclen)
  errs[[idx]] <- z; idx <- idx + 1L
  
  for (i in 1L:veclen) {
    v <- z; v[i] <- 1L; errs[[idx]] <- v; idx <- idx + 1L
  }
  for (i in 1L:(veclen-1L)) for (j in (i+1L):veclen) {
    v <- z; v[i] <- 1L; v[j] <- 1L; errs[[idx]] <- v; idx <- idx + 1L
  }
  for (i in 1L:(veclen-2L)) for (j in (i+1L):(veclen-1L)) for (k in (j+1L):veclen) {
    v <- z; v[i] <- 1L; v[j] <- 1L; v[k] <- 1L; errs[[idx]] <- v; idx <- idx + 1L
  }
  errs
}

## Syndrome lookup table: key = 12-bit syndrome as a string, value = 24-bit error vec
.build_syndrome_lut <- function(H) {
  lut <- new.env(parent = emptyenv())
  for (err in .make_3bit_errors(24L)) {
    syn <- as.integer((H %*% err) %% 2)
    key <- paste(syn, collapse = "")
    lut[[key]] <- err
  }
  lut
}

.DEFAULT_SYNDROME_LUT <- .build_syndrome_lut(DEFAULT_H)

## --- Helpers: sequence<->bits ------------------------------------------------

.seq_to_bits <- function(seq, nt_to_bits = DECODE_GOLAY_NT_TO_BITS) {
  nts <- strsplit(toupper(seq), "", fixed = TRUE)[[1]]
  bits_each <- vapply(nts, function(n) {
    v <- nt_to_bits[[n]]
    if (is.null(v)) stop(sprintf("Invalid nucleotide '%s' in '%s'", n, seq), call. = FALSE)
    v
  }, character(1))
  bit_chars <- strsplit(paste(bits_each, collapse = ""), "", fixed = TRUE)[[1]]
  as.integer(bit_chars)
}

.bits_to_seq <- function(bits, nt_to_bits = ENCODE_GOLAY_NT_TO_BITS) {
  if (length(bits) %% 2L != 0L) stop("Bit vector length must be even.", call. = FALSE)
  bits_to_nt <- setNames(names(nt_to_bits), nt_to_bits)  # e.g., "11" -> "A"
  pairs <- vapply(seq.int(1L, length(bits), by = 2L),
                  function(i) paste0(bits[i], bits[i+1]), character(1))
  paste0(bits_to_nt[pairs], collapse = "")
}

## --- Core: decode a 24-bit received vector via syndrome decoding -------------

.decode_bits <- function(rec_bits, H = DEFAULT_H, lut = .DEFAULT_SYNDROME_LUT) {
  syn <- as.integer((H %*% rec_bits) %% 2)
  key <- paste(syn, collapse = "")
  err <- lut[[key]]
  if (is.null(err)) {
    return(list(corrected = NULL, num_errors = 4L))  # uncorrectable (likely 4-bit)
  }
  corrected <- as.integer((rec_bits + err) %% 2)
  list(corrected = corrected, num_errors = sum(err))
}

## --- Public API: decode Golay 12-nt sequences --------------------------------

#' Decode 12-nt DNA Golay barcodes (vectorized).
#'
#' @param x Character string, character vector, or Biostrings::DNAString/ DNAStringSet.
#' @param nt_to_bits Optional named character vector mapping nts to 2-bit strings.
#'                   Defaults allow N during decoding ("N" -> "11").
#' @return Data frame with columns:
#'         input, corrected (NA if uncorrectable), num_bit_errors (int),
#'         uncorrectable (logical).
correctGolay12 <- function(x, nt_to_bits = DECODE_GOLAY_NT_TO_BITS) {
  # Accept DNAString / DNAStringSet transparently if Biostrings is present
  if (inherits(x, "DNAString") || inherits(x, "DNAStringSet")) x <- as.character(x)
  if (length(x) == 0L) return(data.frame(input=character(), corrected=character(),
                                         num_bit_errors=integer(), uncorrectable=logical(),
                                         stringsAsFactors = FALSE))
  if (!is.character(x)) stop("Input must be character or DNAString(Set).", call. = FALSE)
  
  decode_one <- function(s) {
    if (is.na(s) || nchar(s) == 0L) {
      return(c(corrected = NA_character_, num_bit_errors = NA_integer_, uncorrectable = TRUE))
    }
    if (nchar(s) != 12L) {
      return(c(corrected = NA_character_, num_bit_errors = NA_integer_, uncorrectable = TRUE))
    }
    rec_bits <- .seq_to_bits(s, nt_to_bits = nt_to_bits)
    out <- .decode_bits(rec_bits, H = DEFAULT_H, lut = .DEFAULT_SYNDROME_LUT)
    if (is.null(out$corrected)) {
      c(corrected = NA_character_, num_bit_errors = 4L, uncorrectable = TRUE)
    } else {
      corrected_seq <- .bits_to_seq(out$corrected, nt_to_bits = ENCODE_GOLAY_NT_TO_BITS)
      c(corrected = corrected_seq, num_bit_errors = out$num_errors, uncorrectable = FALSE)
    }
  }
  
  res <- t(vapply(x, decode_one,
                  FUN.VALUE = c(corrected = "", num_bit_errors = 0, uncorrectable = FALSE)))
  res <- as.data.frame(res, stringsAsFactors = FALSE)
  res$num_bit_errors <- as.integer(res$num_bit_errors)
  res$uncorrectable <- as.logical(res$uncorrectable)
  data.frame(input = x, res, row.names = NULL, check.names = FALSE)
}
