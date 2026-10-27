# Subexponential Weighted Counting and Exact Optimisation

Author: Deep Bhattacharjee. Source and verification code for the paper of that
title (48 pages, amsart). Release v1.3.0.

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
  bounded-degree graphs (Theorem 6.15, Corollary 6.16);
- the density of states of Ising and Potts systems with polynomially bounded
  couplings (Theorem 6.18);
- dominating sets, matchings, hitting sets and set splitting (Theorem 6.20),
  sparse permanents over the integers (Theorem 6.21) and counting under a
  congruence (Proposition 6.23);

and, with the exponent n (log log n)^2 / (log n)^{1/4} after a branching
rule,

- the weighted independence polynomial of every graph (Theorem 6.7), hence
  the independence number, the clique number and their counts
  (Corollary 6.8); independent sets under any fixed number of vertex-weight
  budgets on every graph (Theorem 6.10), hence the knapsack problem with an
  arbitrary conflict graph and polynomially bounded weights, and its
  multidimensional form (Corollary 6.11); and Vertex Cover in time
  2^{o(k)} + poly(n) (Theorem 6.13).

Section 7 works out further families, each with a complete proof: several
linear statistics at once with equations and congruences imposed (Theorem
7.1), Max-k-Lin mod p, Max-q-Cut and balanced q-partitions (Corollary 7.4),
f-factors, orientations with prescribed indegrees, proper edge colourings and
nowhere-zero flows on bounded-degree graphs (Theorem 7.5), hafnians of sparse
matrices (Corollary 7.6), the Tutte polynomial on every hyperbola
(x-1)(y-1) = q with q a positive integer, including the chromatic and flow
polynomials at q and the Potts model (Theorem 7.10), the weight enumerator,
minimum distance and maximum-likelihood decoding of codes with sparse parity
checks (Theorem 7.12), matchings, transversals and colourings of
bounded-degree hypergraphs (Theorem 7.14), the output profile of linear-size
circuits (Theorem 7.15), 0-1 and bounded integer programs with O(n) rows of
boundedly many variables (coefficients of any size) and a fixed number of dense
rows (Theorem 7.17), including exact cover with bounded element frequency and
the partially ordered knapsack problem (Corollary 7.18), and growing density:
with parameters that follow the degree ratio and the entropy bound for the
Taylor count (Lemma 3.4), the exponent is O(n g(θ)) with
θ = ck sqrt(log log n / log n) and g(x) = x log^{3/2}(1/x) (Theorems 7.19,
7.20), so Max-Cut, the Ising and Potts ground states, Max-k-SAT and sparse
integer programs at every average degree o(sqrt(log n / log log n)) run in time
2^{o(n)} (Corollaries 7.21, 7.23), and so does the cut profile of every graph
whose complement has that many edges, which includes graphs with Θ(n²) edges
(Corollary 7.22). Remark 7.24 shows that no design for the routine reaches
total degree (n+1) log n / 4.

Section 8 proves k-SAT ∈ SE from the sparsification lemma and SNP ⊆ SE.
Section 9 lists what the method does not give: Max-Cut on general graphs with Θ(n²) edges, the
chromatic and Tutte polynomials as polynomials in q, the quadratic-time
problems of fine-grained complexity, and P versus NP, on which the paper
proves nothing.

## Layout

- `paper/`: LaTeX source (`main.tex` and the section files), figures in
  `paper/figures/`.
- `figures/make_figures.py`: regenerates the four figures (matplotlib).
- `code/check_interpolation.py`: the exact-arithmetic check of Appendix A
  (Python standard library only, `Fraction`; about twenty minutes). Its
  recorded output is `code/check_output.txt`; the expected last line is
  `7 checks passed, 0 failed`.
- `code/check_density.py`: the exact checks of the growing-density section
  (entropy bound, complement transfer on all graphs with at most six vertices,
  the parameter inequalities of Theorem 7.19, Remark 7.24 and the footnote to
  Theorem 7.1; a few seconds). Recorded output `code/check_density_output.txt`,
  expected last line `6 checks passed, 0 failed`.
- `dist/`: the compiled paper, the TeX zip and the arXiv tarball attached to the
  release.
- `.github/workflows/check.yml` runs the check program and compiles the paper;
  `.github/workflows/release.yml` publishes a release by hand.

## Build

    cd paper && latexmk -pdf main.tex
    cd code && python3 check_interpolation.py && python3 check_density.py

## Archive

| Release | DOI |
| --- | --- |
| v1.3.0 | [10.5281/zenodo.23706145](https://doi.org/10.5281/zenodo.23706145) |
| v1.2.0 | [10.5281/zenodo.23698450](https://doi.org/10.5281/zenodo.23698450) |
| v1.1.0 | [10.5281/zenodo.23697614](https://doi.org/10.5281/zenodo.23697614) |
| v1.0.0 | [10.5281/zenodo.23695836](https://doi.org/10.5281/zenodo.23695836) (concept DOI of all versions [10.5281/zenodo.23695835](https://doi.org/10.5281/zenodo.23695835)) |

The code is under the MIT licence (see `LICENSE`); the paper is the author's.
