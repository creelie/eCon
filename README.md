# Subexponential Weighted Counting and Exact Optimisation

Author: Deep Bhattacharjee. Source and verification code for the paper of that
title (31 pages, amsart). Release v1.0.0.

## What the paper proves

OpenAI's preprint *The Exponential Time Hypothesis is False* (25 October 2026)
gives a deterministic algorithm that counts the satisfying assignments of a
k-CNF formula in time 2^{o(n)}. The paper reproduces, with credit and in full,
the part of that algorithm that handles formulas with O(n) clauses (Sections 3
and 4), observes that every step of it is linear over the coefficient ring, and
runs it over truncated polynomial rings and group algebras. From this it proves,
with complete proofs, deterministic algorithms of running time
poly(n, Ω) 2^{O(n (log log n)^2 / sqrt(log n))} for

- the whole weight profile of a constraint system with O(n) constraints of
  bounded scope and integer weights at most Ω (Theorem 6.1), hence Max-k-SAT,
  Min-k-SAT and the number of optimal assignments (Corollary 6.2);
- the number of satisfying assignments of each Hamming weight (Corollary 6.3);
- Max-Cut and the maximum and minimum bisection of a graph with O(n) edges,
  with their counts (Corollary 6.4);
- partition functions, homomorphism counts and discrete graphical models on
  bounded-degree graphs (Theorem 6.13, Corollary 6.14);
- the density of states of Ising and Potts systems with polynomially bounded
  couplings (Theorem 6.16);
- dominating sets, matchings, hitting sets and set splitting (Theorem 6.18),
  sparse permanents over the integers (Theorem 6.19) and counting under a
  congruence (Proposition 6.21);

and, with the exponent n (log log n)^2 / (log n)^{1/4} after a branching
rule,

- the weighted independence polynomial of every graph (Theorem 6.7), hence
  the independence number, the clique number and their counts
  (Corollary 6.8), and Vertex Cover in time 2^{o(k)} + poly(n) (Theorem 6.11).

Section 7 proves k-SAT ∈ SE from the sparsification lemma and SNP ⊆ SE.
Section 8 lists what the method does not give: dense objectives, the
quadratic-time problems of fine-grained complexity, and P versus NP, on which
the paper proves nothing.

## Layout

- `paper/`: LaTeX source (`main.tex` and the section files), figures in
  `paper/figures/`.
- `figures/make_figures.py`: regenerates the three figures (matplotlib).
- `code/check_interpolation.py`: the exact-arithmetic check of Appendix A
  (Python standard library only, `Fraction`; about ten minutes). Its recorded
  output is `code/check_output.txt`; the expected last line is
  `4 checks passed, 0 failed`.
- `dist/`: the compiled paper, the TeX zip and the arXiv tarball attached to the
  release.
- `.github/workflows/check.yml` runs the check program and compiles the paper;
  `.github/workflows/release.yml` publishes a release by hand.

## Build

    cd paper && latexmk -pdf main.tex
    cd code && python3 check_interpolation.py

## Archive

| Release | DOI |
| --- | --- |
| v1.0.0 | [10.5281/zenodo.23695836](https://doi.org/10.5281/zenodo.23695836) (concept DOI of all versions [10.5281/zenodo.23695835](https://doi.org/10.5281/zenodo.23695835)) |

The code is under the MIT licence (see `LICENSE`); the paper is the author's.
