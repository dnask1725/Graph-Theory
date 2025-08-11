# ---------- Synthetic FASTA of amplicon reads  ----------
fasta_txt <- "
>Read_01
ACGTTGACCTGAACT
>Read_02
TGAACTGGTAACTGA
>Read_03
GGTAACTGACGTTGA
>Read_04
ACTGACGTTGACCTG
>Read_05
CCTGAACTGGTAACT
"

# --------- Minimal FASTA parser -------------
read_fasta_text <- function(txt) {
  lines <- strsplit(txt, "\n")[[1]]
  lines <- trimws(lines)
  lines <- lines[nchar(lines) > 0]
  headers <- grep("^>", lines)
  seq_names <- sub("^>", "", lines[headers])
  seq_values <- character(length(seq_names))
  for (i in seq_along(headers)) {
    start <- headers[i] + 1
    end   <- if (i < length(headers)) headers[i+1] - 1 else length(lines)
    if (start <= end) seq_values[i] <- paste0(lines[start:end], collapse = "")
  }
  names(seq_values) <- seq_names
  seq_values
}

# --------- Build overlap graph edges: suffix_k(s) == prefix_k(t) -------
overlap_edges <- function(seq_vec, k = 3) {
  nm <- names(seq_vec)
  A  <- as.character(seq_vec)
  # Precompute prefixes/suffixes
  pref <- substr(A, 1, k)
  suff <- substr(A, pmax(1, nchar(A) - k + 1), nchar(A))
  edges_from <- character(0)
  edges_to   <- character(0)
  for (i in seq_along(A)) {
    # match suffix of i to prefix of everyone else
    hits <- which(pref == suff[i] & seq_along(A) != i)
    if (length(hits)) {
      edges_from <- c(edges_from, rep(nm[i], length(hits)))
      edges_to   <- c(edges_to,   nm[hits])
    }
  }
  data.frame(from = edges_from, to = edges_to, stringsAsFactors = FALSE)
}

# ---------------- Run on the synthetic reads (k = 3) -------------------
reads <- read_fasta_text(fasta_txt)
edges <- overlap_edges(reads, k = 3)
edges
