# Needleman–Wunsch Global Sequence Alignment (R)

This repository contains an implementation of the **Needleman–Wunsch algorithm**
for global sequence alignment written in **R**.

The algorithm uses dynamic programming to compute optimal global alignments
between two nucleotide sequences and supports **multiple optimal alignments**
via recursive traceback.

---

## 🧬 Features
- Global sequence alignment (Needleman–Wunsch)
- Customizable scoring scheme (match, mismatch, gap)
- Supports multiple optimal alignments
- Outputs alignment score and full scoring matrix
- Written for clarity and educational use

---

## 🚀 Usage

```r
source("needleman_wunsch_global_alignment.R")

result <- global_alignment(
  seq1 = "ATCG",
  seq2 = "ATG",
  match = 1,
  mismatch = -1,
  gap = -2
)

result$score
result$alignments
result$score_matrix

##Example Output
Sequence 1: ATCG
Sequence 2: ATG

Optimal Alignment Score: 0
Number of optimal alignments: 1

Seq1: ATCG
Seq2: AT-G
      || |

## ⏱️ Runtime & Space Complexity

Let **n** be the length of sequence 1 and **m** be the length of sequence 2.

### Time Complexity
- Filling the dynamic programming matrix requires evaluating three scores
  (diagonal, up, left) for each cell.
- This results in a time complexity of **O(n × m)**.

- The traceback step may explore multiple paths when ties occur in the scoring
  matrix.
  - In the **worst case**, the number of optimal alignments can grow
    exponentially.
  - However, for most biological sequences and scoring schemes, the traceback
    runtime is much smaller than matrix construction.

### Space Complexity
- The scoring matrix requires **O(n × m)** space.
- The traceback matrix also requires **O(n × m)** space.
- Overall space complexity is **O(n × m)**.

### Practical Notes
- Global alignment guarantees an optimal alignment across the full lengths of
  both sequences, but at the cost of quadratic time and space.
- For very long sequences (e.g., whole genomes), heuristic or local alignment
  methods (such as BLAST or Smith–Waterman) are often preferred.
