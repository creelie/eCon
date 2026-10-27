"""Exact checks for the growing-density part of the paper (Section 7.9).

Five groups, all in integer or rational arithmetic:
  (A) the entropy bound of Lemma 3.4 in the form C_q q^q (n-q)^(n-q) <= 2^q n^n;
  (B) the complement identity (7.7) behind Corollary 7.22, on every graph with
      at most six vertices and on weighted instances;
  (C) the parameter inequalities in the proof of Theorem 7.19, at sample values
      of lambda and of the degree ratio a, in the logarithmic domain;
  (D) the top coefficient (-2)^m (W-1)^m of the cut-factor product in
      Remark 7.24;
  (E) the two values in the footnote to the proof of Theorem 7.1.
The asymptotic statements rest on their written proofs, not on these checks.
"""
from fractions import Fraction
from itertools import combinations
from math import comb, isqrt
import random
import sys

passed = 0
failed = 0


def report(name, ok):
    global passed, failed
    if ok:
        passed += 1
    else:
        failed += 1
    print(("ok    " if ok else "FAIL  ") + name)


def ceil_log2(x):
    return (x - 1).bit_length()


# (A) Lemma 3.4
def check_entropy():
    ok = True
    count = 0
    for n in range(1, 201):
        for q in range(0, n // 2 + 1):
            cq = sum((1 << s) * comb(n, s) for s in range(q + 1))
            lhs = cq * q ** q * (n - q) ** (n - q)
            ok &= lhs <= (1 << q) * n ** n
            count += 1
    report("(A) entropy bound for 1 <= n <= 200, 0 <= q <= n/2 (%d cases)" % count, ok)


# (B) complement identity and profile transfer
def cut_profile(n, weight):
    prof = {}
    for mask in range(1 << n):
        size = bin(mask).count("1")
        cut = 0
        for (u, v), w in weight.items():
            if ((mask >> u) & 1) != ((mask >> v) & 1):
                cut += w
        prof[(size, cut)] = prof.get((size, cut), 0) + 1
    return prof


def check_complement():
    ok = True
    graphs = 0
    for n in range(0, 7):
        pairs = list(combinations(range(n), 2))
        for gmask in range(1 << len(pairs)):
            g = {e: 1 for j, e in enumerate(pairs) if (gmask >> j) & 1}
            h = {e: 1 for j, e in enumerate(pairs) if not (gmask >> j) & 1}
            pg = cut_profile(n, g)
            ph = cut_profile(n, h)
            moved = {}
            for (i, j), c in ph.items():
                moved[(i, i * (n - i) - j)] = moved.get((i, i * (n - i) - j), 0) + c
            ok &= moved == pg
            graphs += 1
    report("(B1) cut profile of every graph on at most 6 vertices from its complement (%d graphs)" % graphs, ok)
    rng = random.Random(20261027)
    ok = True
    for trial in range(200):
        n = rng.randint(2, 7)
        omega0 = rng.randint(1, 5)
        pairs = list(combinations(range(n), 2))
        hedges = [e for e in pairs if rng.random() < 0.3]
        omega = {e: omega0 for e in pairs}
        eta = {}
        for e in hedges:
            omega[e] = rng.randint(0, omega0)
            eta[e] = omega0 - omega[e]
        pw = cut_profile(n, omega)
        pe = cut_profile(n, eta)
        moved = {}
        for (i, j), c in pe.items():
            key = (i, omega0 * i * (n - i) - j)
            moved[key] = moved.get(key, 0) + c
        ok &= moved == pw
    report("(B2) weighted form of the identity on 200 random instances", ok)


# (C) parameters of Theorem 7.19
def params(lam, a):
    nu = ceil_log2(lam)
    t = 0
    while 4 ** (t + 1) * a * a * nu <= lam:
        t += 1
    r = 1 << t
    Q = lam // (128 * t)
    R = 32 + isqrt(256 * Q * nu - 1) + 1  # 32 + ceil(16 sqrt(Q nu))
    return nu, t, r, Q, R


def check_params():
    ok = True
    cases = 0
    for e in (40, 48, 64, 80, 100, 128):
        lam = 1 << e
        nu = ceil_log2(lam)
        amax = isqrt(lam // (nu << 30))  # theta^2 = a^2 nu / lam <= 2^-30
        for a in sorted({1, 2, 7, max(1, amax // 1000), max(1, amax // 3), amax}):
            if a < 1 or a * a * nu * (1 << 30) > lam:
                continue
            nu, t, r, Q, R = params(lam, a)
            c = []
            c.append(4 ** t * a * a * nu <= lam < 4 ** (t + 1) * a * a * nu)  # t = floor(L)
            c.append(t >= 15 and 2 * t <= nu)
            c.append(lam >= 256 * t and 2 * Q * 128 * t >= lam)  # Q >= lam/(256 t)
            c.append((R - 32 - 1) ** 2 < 256 * Q * nu <= (R - 32) ** 2)  # R is the ceiling
            c.append(R >= 32 and 6 * R <= Q)
            logD = 2 * Q * t
            c.append(64 * logD <= lam)
            # 2(b+w)D^2 <= 4 lam^2 2^(lam/32) <= 2^(lam-1), with b, w <= lam^2
            c.append(ceil_log2(4 * lam * lam) + 2 * logD <= lam - 1)
            c.append(ceil_log2(2 * Q) + 2 * t <= lam - 1)  # 2 Q r^2 < n
            c.append(R * R >= 256 * Q * nu)  # Hoeffding exponent >= 32 nu
            c.append(r * r <= 1 << nu)
            # chi = aR/Q <= 1/10 and, squared, chi^2 <= 784^2 * theta^2 * t
            c.append(10 * a * R <= Q)
            c.append((a * R) ** 2 * lam <= 784 ** 2 * a * a * nu * t * Q * Q)
            ok &= all(c)
            cases += 1
    report("(C) parameter inequalities of Theorem 7.19 at %d sample points" % cases, ok)


# (D) top coefficient of the cut-factor product
def poly_mul(p, q):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = tuple(sorted(k1 + k2))
            out[k] = out.get(k, 0) + c1 * c2
    return {k: c for k, c in out.items() if c != 0}


def check_topcoeff():
    ok = True
    for m in range(1, 6):
        for W in (Fraction(2), Fraction(3), Fraction(-1, 2)):
            prod = {(): Fraction(1)}
            for e in range(m):
                u, v = ("u%d" % e,), ("v%d" % e,)
                factor = {(): Fraction(1), u: W - 1, v: W - 1, tuple(sorted(u + v)): -2 * (W - 1)}
                prod = poly_mul(prod, factor)
            top = max(prod, key=len)
            ok &= len(top) == 2 * m and prod[top] == (-2) ** m * (W - 1) ** m
    report("(D) degree 2m and top coefficient (-2)^m (W-1)^m for m <= 5", ok)


# (E) the footnote values
def check_footnote():
    W1, W2, Y = Fraction(2), Fraction(3), Fraction(1, 2)
    product = (1 + (W1 - 1) * Y) * (1 + (W2 - 1) * Y)
    psi = 1 + (W1 * W2 - 1) * Y
    report("(E) product 3 and Psi_v 7/2 at W1=2, W2=3, Y=1/2", product == 3 and psi == Fraction(7, 2))


if __name__ == "__main__":
    check_entropy()
    check_complement()
    check_params()
    check_topcoeff()
    check_footnote()
    print("%d checks passed, %d failed" % (passed, failed))
    sys.exit(1 if failed else 0)
