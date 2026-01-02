############################################################
# Needleman–Wunsch Global Sequence Alignment
# Author: Jordan Scollick-Karnitz
# Language: R
#
# Description:
# Implements the Needleman–Wunsch algorithm for global
# sequence alignment using dynamic programming.
# Supports multiple optimal alignments via traceback.
############################################################

global_alignment <- function(seq1 = "ATCG",
                             seq2 = "ATG",
                             match = 1,
                             mismatch = -1,
                             gap = -2) {
  
  # Sequence lengths
  n <- nchar(seq1)
  m <- nchar(seq2)
  
  # Split sequences into individual characters
  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]
  
  # Initialize score matrix (n+1 x m+1)
  F <- matrix(0, nrow = n + 1, ncol = m + 1)
  
  # Initialize first column with gap penalties
  for (i in 2:(n + 1)) {
    F[i, 1] <- F[i - 1, 1] + gap
  }
  
  # Initialize first row with gap penalties
  for (j in 2:(m + 1)) {
    F[1, j] <- F[1, j - 1] + gap
  }
  
  # Traceback matrix to store optimal directions
  traceback <- vector("list", (n + 1) * (m + 1))
  dim(traceback) <- c(n + 1, m + 1)
  
  # Fill scoring and traceback matrices
  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      
      score_diag <- F[i - 1, j - 1] +
        ifelse(s1[i - 1] == s2[j - 1], match, mismatch)
      
      score_up   <- F[i - 1, j] + gap
      score_left <- F[i, j - 1] + gap
      
      best_score <- max(score_diag, score_up, score_left)
      F[i, j] <- best_score
      
      directions <- character(0)
      if (score_diag == best_score) directions <- c(directions, "diag")
      if (score_up   == best_score) directions <- c(directions, "up")
      if (score_left == best_score) directions <- c(directions, "left")
      
      traceback[[i, j]] <- directions
    }
  }
  
  final_score <- F[n + 1, m + 1]
  alignments <- list()
  
  # Recursive traceback to recover all optimal alignments
  trace <- function(i, j, aln1 = "", aln2 = "") {
    
    if (i == 1 && j == 1) {
      alignments <<- append(alignments, list(c(aln1, aln2)))
      return()
    }
    
    if (i == 1) {
      trace(i, j - 1,
            paste0("-", aln1),
            paste0(s2[j - 1], aln2))
      return()
    }
    
    if (j == 1) {
      trace(i - 1, j,
            paste0(s1[i - 1], aln1),
            paste0("-", aln2))
      return()
    }
    
    for (d in traceback[[i, j]]) {
      if (d == "diag") {
        trace(i - 1, j - 1,
              paste0(s1[i - 1], aln1),
              paste0(s2[j - 1], aln2))
      } else if (d == "up") {
        trace(i - 1, j,
              paste0(s1[i - 1], aln1),
              paste0("-", aln2))
      } else if (d == "left") {
        trace(i, j - 1,
              paste0("-", aln1),
              paste0(s2[j - 1], aln2))
      }
    }
  }
  
  trace(n + 1, m + 1)
  
  return(list(
    score = final_score,
    alignments = alignments,
    score_matrix = F
  ))
}

############################################################
# Example usage
############################################################

cat("GLOBAL SEQUENCE ALIGNMENT (Needleman–Wunsch)\n\n")

result <- global_alignment("ATCG", "ATG")

cat("Sequence 1: ATCG\n")
cat("Sequence 2: ATG\n\n")
cat("Scoring Scheme:\n")
cat("  Match: +1\n")
cat("  Mismatch: -1\n")
cat("  Gap: -2\n\n")
cat("Optimal Alignment Score:", result$score, "\n\n")
cat("Number of optimal alignments:", length(result$alignments), "\n\n")

for (i in seq_along(result$alignments)) {
  aln <- result$alignments[[i]]
  cat("Alignment", i, ":\n")
  cat("Seq1:", aln[1], "\n")
  cat("Seq2:", aln[2], "\n")
  
  matches <- ifelse(
    strsplit(aln[1], "")[[1]] ==
      strsplit(aln[2], "")[[1]],
    "|", " "
  )
  cat("      ", paste(matches, collapse = ""), "\n\n")
}

cat("Score Matrix:\n")
print(result$score_matrix)
