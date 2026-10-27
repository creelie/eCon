#!/usr/bin/env python3
"""Exact-arithmetic check of the shared-summary interpolation algorithm
(weighted form) on small instances.

The algorithm is the one of OpenAI's preprint "The Exponential Time
Hypothesis is False" (October 25, 2026), Sections 3-4, run here with
coefficients in the truncated polynomial ring A = Q[W]/(W^N) so that it
returns generating polynomials (Max-SAT profile, independence polynomial,
weight enumerator of a code) and not only counts.  Every number is a Python Fraction; nothing is
floating point.  The design parameters are tiny (r = 2, Q = 3, R = 2,
D = 64), far from the asymptotic regime, because correctness of the
identity H(0) = Z needs only Definition 3.1, not the size bounds.
"""
from fractions import Fraction as Fr
from itertools import product, combinations
import random, sys

# ---------------------------------------------------------------- algebra
class Poly:
    """Element of Q[W]/(W^N) (coefficient list of length N)."""
    N = 1
    __slots__ = ("c",)
    def __init__(self, c=None):
        if c is None:
            c = [Fr(0)] * Poly.N
        self.c = c
    @staticmethod
    def const(a):
        c = [Fr(0)] * Poly.N; c[0] = Fr(a); return Poly(c)
    @staticmethod
    def W():
        c = [Fr(0)] * Poly.N
        if Poly.N > 1: c[1] = Fr(1)
        return Poly(c)
    def __add__(s, o):
        if not isinstance(o, Poly): o = Poly.const(o)
        return Poly([a + b for a, b in zip(s.c, o.c)])
    __radd__ = __add__
    def __neg__(s): return Poly([-a for a in s.c])
    def __sub__(s, o):
        if not isinstance(o, Poly): o = Poly.const(o)
        return s + (-o)
    def __rsub__(s, o): return Poly.const(o) - s
    def __mul__(s, o):
        if not isinstance(o, Poly):
            o = Fr(o); return Poly([a * o for a in s.c])
        N = Poly.N; out = [Fr(0)] * N
        for i, a in enumerate(s.c):
            if a == 0: continue
            for j in range(N - i):
                b = o.c[j]
                if b: out[i + j] += a * b
        return Poly(out)
    __rmul__ = __mul__
    def iszero(s): return all(a == 0 for a in s.c)
    def __eq__(s, o): return (s - o).iszero()
    def __repr__(s):
        return " + ".join(f"{a}*W^{i}" for i, a in enumerate(s.c) if a) or "0"

class ZPoly:
    """Element of A[z]/(z^J) with A = Poly; list of J Poly's."""
    J = 1
    __slots__ = ("c",)
    def __init__(s, c=None):
        s.c = c if c is not None else [Poly() for _ in range(ZPoly.J)]
    @staticmethod
    def lift(a):  # a in A or scalar
        if not isinstance(a, Poly): a = Poly.const(a)
        c = [Poly() for _ in range(ZPoly.J)]; c[0] = a; return ZPoly(c)
    @staticmethod
    def from_qz(coeffs):  # list of Fractions in z
        c = [Poly() for _ in range(ZPoly.J)]
        for j, a in enumerate(coeffs[:ZPoly.J]): c[j] = Poly.const(a)
        return ZPoly(c)
    def __add__(s, o):
        if not isinstance(o, ZPoly): o = ZPoly.lift(o)
        return ZPoly([a + b for a, b in zip(s.c, o.c)])
    __radd__ = __add__
    def __neg__(s): return ZPoly([-a for a in s.c])
    def __sub__(s, o):
        if not isinstance(o, ZPoly): o = ZPoly.lift(o)
        return s + (-o)
    def __rsub__(s, o): return ZPoly.lift(o) - s
    def __mul__(s, o):
        if not isinstance(o, ZPoly): o = ZPoly.lift(o)
        J = ZPoly.J; out = [Poly() for _ in range(J)]
        for i, a in enumerate(s.c):
            if a.iszero(): continue
            for j in range(J - i):
                b = o.c[j]
                if not b.iszero(): out[i + j] = out[i + j] + a * b
        return ZPoly(out)
    __rmul__ = __mul__

# ------------------------------------------------------------ rational poly in T
def padd(p, q):
    n = max(len(p), len(q)); p = p + [Fr(0)] * (n - len(p)); q = q + [Fr(0)] * (n - len(q))
    return [a + b for a, b in zip(p, q)]
def pmul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): out[i + j] += a * b
    return out
def pscale(p, s): return [a * s for a in p]
def peval(p, t): return sum((a * Fr(t) ** i for i, a in enumerate(p)), Fr(0))
def pshift(p, t):
    """coefficients of p(t+z) in z"""
    out = [Fr(0)] * len(p)
    for i, a in enumerate(p):
        if a == 0: continue
        for j in range(i + 1):
            out[j] += a * binom(i, j) * Fr(t) ** (i - j)
    return out
def binom(n, k):
    from math import comb
    return Fr(comb(n, k))

# ---------------------------------------------------------------- linear algebra
def rank(rows):
    rows = [r[:] for r in rows]; rk = 0; ncol = len(rows[0]) if rows else 0
    for col in range(ncol):
        piv = next((i for i in range(rk, len(rows)) if rows[i][col] != 0), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        pv = rows[rk][col]
        rows[rk] = [a / pv for a in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][col] != 0:
                f = rows[i][col]; rows[i] = [a - f * b for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk

def solve(Mat, rhs):
    """Solve Mat h = rhs over Q with rhs entries in A (Poly). Mat square invertible."""
    n = len(Mat); M = [r[:] for r in Mat]; b = rhs[:]
    for col in range(n):
        piv = next(i for i in range(col, n) if M[i][col] != 0)
        M[col], M[piv] = M[piv], M[col]; b[col], b[piv] = b[piv], b[col]
        pv = M[col][col]; M[col] = [a / pv for a in M[col]]; b[col] = b[col] * (1 / pv)
        for i in range(n):
            if i != col and M[i][col] != 0:
                f = M[i][col]; M[i] = [a - f * c for a, c in zip(M[i], M[col])]; b[i] = b[i] - b[col] * f
    return b

# ---------------------------------------------------------------- the design
class Design:
    def __init__(self, n, scopes, r=2, Q=3, R=2, A0=0, seed=1):
        assert R < Q
        self.n, self.scopes, self.r, self.Q, self.R = n, scopes, r, Q, R
        rng = random.Random(seed)
        seqs = list(product(range(1, r + 1), repeat=Q))          # [r]^Q
        self.B = [(a, d) for a in seqs for d in seqs]              # labels
        D = len(self.B); self.D = D; idx = {u: i for i, u in enumerate(self.B)}
        # swap permutations Pi_s as index maps
        def theta(u, s):
            a, d = list(u[0]), list(u[1]); a[s], d[s] = d[s], a[s]; return (tuple(a), tuple(d))
        self.Pi = [[idx[theta(u, s)] for u in self.B] for s in range(Q)]   # (Pi_s zeta)_u = zeta_{theta_s(u)}
        # eigenbasis
        funcs = []  # (function on [r]^2 as dict, eigenvalue)
        for al in range(1, r + 1): funcs.append(({(al, al): Fr(1)}, +1))
        for al in range(1, r + 1):
            for de in range(al + 1, r + 1):
                funcs.append(({(al, de): Fr(1), (de, al): Fr(1)}, +1))
                funcs.append(({(al, de): Fr(1), (de, al): Fr(-1)}, -1))
        assert len(funcs) == r * r
        I = set(range(A0, A0 + R)); self.I = I
        basis = []  # (vector over B, j)
        for choice in product(range(len(funcs)), repeat=Q):
            j = sum(1 for c in choice if funcs[c][1] == -1)
            vec = [Fr(1)] * D
            for i, u in enumerate(self.B):
                for s in range(Q):
                    vec[i] *= funcs[choice[s]][0].get((u[0][s], u[1][s]), Fr(0))
            basis.append((vec, j))
        assert rank([v for v, _ in basis]) == D
        self.CB = [v for v, j in basis if j not in I]      # measurement rows
        # W0 = 1/2 sum (I - Pi_s), V0 = prod Pi_s  as matrices (D x D)
        def permmat(pi):
            M = [[Fr(0)] * D for _ in range(D)]
            for u in range(D): M[u][pi[u]] = Fr(1)
            return M
        def matmul(X, Y):
            return [[sum(X[i][k] * Y[k][j] for k in range(D) if X[i][k]) for j in range(D)] for i in range(D)]
        Id = [[Fr(int(i == j)) for j in range(D)] for i in range(D)]
        W0 = [[Fr(0)] * D for _ in range(D)]
        for s in range(Q):
            P = permmat(self.Pi[s])
            for i in range(D):
                for j in range(D): W0[i][j] += (Id[i][j] - P[i][j]) / 2
        V0 = Id
        for s in range(Q): V0 = matmul(V0, permmat(self.Pi[s]))
        # f_I(W0) = (-1)^A0 sum_{j<R} (-2)^j binom(W0 - A0, j)
        fW = [[Fr(0)] * D for _ in range(D)]
        term = Id
        for j in range(R):
            coef = Fr((-1) ** A0) * Fr(-2) ** j
            for i in range(D):
                for k in range(D): fW[i][k] += coef * term[i][k]
            # term <- term * (W0 - (A0 + j) I) / (j+1)
            Sh = [[W0[i][k] - (Fr(A0 + j) if i == k else 0) for k in range(D)] for i in range(D)]
            term = matmul(term, Sh); term = [[a / (j + 1) for a in row] for row in term]
        self.M0 = matmul(V0, fW)
        # checks: M0 fixes eigenvectors with j in I; (4.4)
        for vec, j in basis:
            img = [sum(self.M0[i][k] * vec[k] for k in range(D)) for i in range(D)]
            if j in I: assert img == vec, "M0 must fix basis vectors with j in I"
        for u in range(D):
            for y in range(D):
                if self.M0[u][y] != 0:
                    mism = sum(1 for s in range(Q) if self.B[y][0][s] != self.B[u][1][s])
                    assert mism <= R - 1, "(4.4) violated"
        # rows P_u(T): dict y -> poly in T
        self.Prow = []
        for u in range(D):
            row = {}
            for y in range(D):
                row[y] = padd(row.get(y, [Fr(0)]), [Fr(int(u == y)) - self.M0[u][y]])
                if self.M0[u][y] != 0:
                    phi = [Fr(1)]
                    for s in range(Q):
                        if self.B[y][0][s] != self.B[u][1][s]:
                            phi = pmul(phi, [Fr(1), Fr(-1, s + 1)])   # (1 - T/(s+1)), points are 1..Q
                    row[y] = padd(row.get(y, [Fr(0)]), pscale(phi, self.M0[u][y]))
            self.Prow.append(row)
        for u in range(D):   # P_u(0) = delta_u
            for y, p in self.Prow[u].items():
                assert peval(p, 0) == (1 if y == u else 0)
        # labels a_v, d_f (random here; the paper derandomises their choice)
        self.a = [rng.choice(seqs) for _ in range(n)]
        self.d = [rng.choice(seqs) for _ in scopes]
        self.A = [(f, v) for f, I_f in enumerate(scopes) for v in I_f]
        label = {(f, v): idx[(self.a[v], self.d[f])] for f, v in self.A}
        self.label = label
        cnt = {}; self.copy = {}
        for fv in self.A:
            u = label[fv]; self.copy[fv] = cnt.get(u, 0); cnt[u] = cnt.get(u, 0) + 1
        self.h = max(cnt.values())
        # P_{f,v}(T) as row over [n]: pos (c,u') occupied by (f',v') with same copy c
        occupied = {(self.copy[fv], label[fv]): fv[1] for fv in self.A}
        self.Pfv = {}
        for fv in self.A:
            c, u = self.copy[fv], label[fv]
            row = [[Fr(0)] for _ in range(n)]
            for y, p in self.Prow[u].items():
                if (c, y) in occupied:
                    vp = occupied[(c, y)]; row[vp] = padd(row[vp], p)
            self.Pfv[fv] = row
            assert [peval(p, 0) for p in row] == [Fr(int(i == fv[1])) for i in range(n)]
        # measurement matrix Chat = diag(CB,...,CB) iota  : (h * cB) x n
        self.Chat = []
        for c in range(self.h):
            for crow in self.CB:
                row = [Fr(0)] * n
                for (cc, y), vp in occupied.items():
                    if cc == c: row[vp] += crow[y]
                self.Chat.append(row)
        self.S = [[[v for v in range(n) if self.a[v][t] == g] for g in range(1, r + 1)] for t in range(Q)]
        self.sigma = [[self.d[f][t] - 1 for t in range(Q)] for f in range(len(scopes))]  # group index 0..r-1
        # verify Definition 3.1 (4): row of P_{f,v}(t) on S_{t,g} lies in rowspace of C_{t,g}
        for t in range(Q):
            for g in range(r):
                S = self.S[t][g]
                if not S: continue
                Ctg = [[row[v] for v in S] for row in self.Chat]
                rk = rank(Ctg) if Ctg else 0
                for fv in self.A:
                    if self.sigma[fv[0]][t] == g: continue
                    prow = [peval(self.Pfv[fv][v], t + 1) for v in S]
                    assert rank(Ctg + [prow]) == rk, "kernel condition fails"
        self.cmeas = len(self.Chat)

    def meas(self, t, g, xS):
        """measurement vector C_{t,g} x|S_{t,g}  (xS dict v->value)"""
        return tuple(sum(row[v] * xS[v] for v in self.S[t][g]) for row in self.Chat)

# ---------------------------------------------------------------- the algorithm
def count_weighted(n, scopes, factors, design, dH):
    """factors[f](inputs: list of ZPoly, one per v in scopes[f]) -> ZPoly.
    Returns H(0) = sum_x prod_f Phi_f(x) in A, via Taylor data at t=1..Q."""
    Q, r = design.Q, design.r
    q = dH // Q; J = q + 1; ZPoly.J = J
    data = {}   # (t, j) -> Poly
    for t in range(Q):
        tt = t + 1
        # P_{f,v}(t) numeric rows and Delta rows (coefficients in z)
        Pt = {fv: [peval(p, tt) for p in row] for fv, row in design.Pfv.items()}
        Dz = {fv: [pshift(p, tt) for p in row] for fv, row in design.Pfv.items()}  # p(t+z) coeffs; Delta = minus constant
        Ht = ZPoly()
        groups = [[f for f in range(len(scopes)) if design.sigma[f][t] == g] for g in range(r)]
        for esize in range(0, min(q, n) + 1):
            for E in combinations(range(n), esize):
                for usize in range(esize + 1):
                    for U in combinations(E, usize):
                        sign = (-1) ** (esize - usize)
                        Uset = set(U)
                        # fixed corrections <Delta P_{f,v}(z), 1_U> in Q[z]/(z^J)
                        corr = {}
                        for fv in design.A:
                            cz = [Fr(0)] * J
                            for v in Uset:
                                coeffs = Dz[fv][v]
                                for j in range(1, min(J, len(coeffs))):
                                    cz[j] += coeffs[j]
                            corr[fv] = cz
                        def Ag(g, x):   # x: full rational vector (list of Fr), returns ZPoly
                            out = ZPoly.lift(1)
                            for f in groups[g]:
                                inputs = []
                                for v in scopes[f]:
                                    base = sum(Pt[(f, v)][w] * x[w] for w in range(n))
                                    cz = corr[(f, v)][:]; cz[0] += base
                                    inputs.append(ZPoly.from_qz(cz))
                                out = out * factors[f](inputs)
                            return out
                        # Lemma 3.2 over attained measurement tuples
                        allx = list(product([0, 1], repeat=n))
                        tuples = {}
                        for x in allx:
                            key = tuple(design.meas(t, g, {v: Fr(x[v]) for v in design.S[t][g]}) for g in range(r))
                            tuples.setdefault(key, x)
                        total = ZPoly()
                        for key, xbar in tuples.items():
                            prod_ = ZPoly.lift(1); empty = False
                            for g in range(r):
                                S = design.S[t][g]
                                sg = ZPoly()
                                for bits in product([0, 1], repeat=len(S)):
                                    xS = {v: Fr(b) for v, b in zip(S, bits)}
                                    if any(xS[v] != 1 for v in E if v in xS): continue
                                    if design.meas(t, g, xS) != key[g]: continue
                                    x = [Fr(xbar[w]) for w in range(n)]
                                    for v in S: x[v] = xS[v]
                                    sg = sg + Ag(g, x)
                                prod_ = prod_ * sg
                            total = total + prod_
                        Ht = Ht + total * Fr(sign)
        for j in range(J): data[(tt, j)] = Ht.c[j]
    # Hermite interpolation: H of degree < QJ
    m = Q * J
    rows, rhs = [], []
    for t in range(1, Q + 1):
        for j in range(J):
            rows.append([binom(i, j) * Fr(t) ** (i - j) if i >= j else Fr(0) for i in range(m)])
            rhs.append(data[(t, j)])
    h = solve(rows, rhs)
    return h[0], h

# ---------------------------------------------------------------- instances
def lit(y, pos):  # literal polynomial: y or 1-y
    return y if pos else (1 - y)

def clause_factor(signs, weighted):
    def phi(inputs):
        prod_ = ZPoly.lift(1)
        for y, s in zip(inputs, signs): prod_ = prod_ * (1 - lit(y, s))
        sat = 1 - prod_
        if not weighted: return sat
        return 1 + sat * (Poly.W() - 1)      # W^{[sat]}
    return phi

def brute_weighted(n, scopes, factors_bool):
    tot = Poly()
    for x in product([0, 1], repeat=n):
        p = Poly.const(1)
        for f, I in enumerate(scopes): p = p * factors_bool[f]([x[v] for v in I])
        tot = tot + p
    return tot

def run_case(name, n, scopes, factors, factors_bool, N, dH, seed):
    Poly.N = N
    des = Design(n, scopes, seed=seed)
    val, _ = count_weighted(n, scopes, factors, des, dH)
    ref = brute_weighted(n, scopes, factors_bool)
    ok = (val == ref)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: n={n} m={len(scopes)} h={des.h} c={des.cmeas} "
          f"Q={des.Q} R={des.R} D={des.D} q={dH // des.Q}")
    print(f"       algorithm: {val}")
    print(f"       brute    : {ref}")
    return ok

def main():
    rng = random.Random(7); results = []
    # case 1: #2-SAT count, A = Q
    n = 6; m = 7
    cls = [(tuple(rng.sample(range(n), 2)), (rng.random() < 0.5, rng.random() < 0.5)) for _ in range(m)]
    scopes = [c[0] for c in cls]
    facs = [clause_factor(c[1], False) for c in cls]
    facs_b = [(lambda s: (lambda xs: Poly.const(1 - int(all((x if not sg else 1 - x) for x, sg in zip(xs, s))))))(c[1]) for c in cls]
    results.append(run_case("#2-SAT count over Q", n, scopes, facs, facs_b, 1, 2 * m * 2, 1))
    # case 2: Max-2-SAT profile sum_x W^{#satisfied clauses}
    facs = [clause_factor(c[1], True) for c in cls]
    def mk(s):
        def fb(xs):
            sat = 1 - int(all((x if not sg else 1 - x) for x, sg in zip(xs, s)))
            c = [Fr(0)] * Poly.N; c[sat] = Fr(1); return Poly(c)
        return fb
    facs_b = [mk(c[1]) for c in cls]
    results.append(run_case("Max-2-SAT profile over Q[W]/(W^(m+1))", n, scopes, facs, facs_b, m + 1, 2 * m * 2, 2))
    # case 3: independence polynomial of a 6-vertex graph (edge factors 1 - y_u y_v, vertex factors 1+(W-1)y_v)
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5), (1, 4)]
    scopes = [list(e) for e in edges] + [[v] for v in range(n)]
    facs = [(lambda ins: 1 - ins[0] * ins[1]) for _ in edges] + [(lambda ins: 1 + ins[0] * (Poly.W() - 1)) for _ in range(n)]
    def eb(xs): return Poly.const(1 - xs[0] * xs[1])
    def vb(xs):
        c = [Fr(0)] * Poly.N; c[xs[0]] = Fr(1); return Poly(c)
    facs_b = [eb] * len(edges) + [vb] * n
    dH = (2 * len(edges) + n) * 2
    results.append(run_case("independence polynomial of C6 plus chord", n, scopes, facs, facs_b, n + 1, dH, 3))
    # case 4: Max-Cut profile sum_x W^{cut(x)} on the same graph
    scopes = [list(e) for e in edges]
    facs = [(lambda ins: 1 + (ins[0] + ins[1] - 2 * ins[0] * ins[1]) * (Poly.W() - 1)) for _ in edges]
    def cb(xs):
        c = [Fr(0)] * Poly.N; c[int(xs[0] != xs[1])] = Fr(1); return Poly(c)
    results.append(run_case("Max-Cut profile of C6 plus chord", n, scopes, facs, [cb] * len(edges), len(edges) + 1, 2 * len(edges) * 2, 4))
    # case 5: weight enumerator of the binary code of length 6 with parity checks
    # {0,1,2}, {2,3,4}, {4,5,0} (Theorem 7.14): factor (1 + prod(1-2y))/2 per check,
    # 1 + (W-1) y_v per bit, over Q[W]/(W^7)
    checks = [[0, 1, 2], [2, 3, 4], [4, 5, 0]]
    scopes = [list(c) for c in checks] + [[v] for v in range(n)]
    def parity_even(ins):
        p = ZPoly.lift(1)
        for y in ins: p = p * (1 - 2 * y)
        return (1 + p) * Fr(1, 2)
    facs = [parity_even for _ in checks] + [(lambda ins: 1 + ins[0] * (Poly.W() - 1)) for _ in range(n)]
    def pb(xs): return Poly.const(1 - (sum(xs) % 2))
    facs_b = [pb] * len(checks) + [vb] * n
    results.append(run_case("weight enumerator of a [6,3] code over Q[W]/(W^7)", n, scopes, facs, facs_b, n + 1, (3 * 3 + n) * 2, 5))
    # case 6: nowhere-zero Z_3-flows of K_{2,3} (Lemma 7.10): edge e carries the value
    # 1 + x_e in {1,2}, oriented from side {a0,a1} to side {b0,b1,b2}; a vertex factor is the
    # multilinear extension of the indicator that the signed sum of its edge values is 0 mod 3
    A, B = [0, 1], [2, 3, 4]
    edges6 = [(a, b) for a in A for b in B]           # 6 edges, variable x_e per edge
    def mle(size, func):
        """multilinear extension of func: {0,1}^size -> Q, as a factor on ZPoly inputs"""
        table = {U: Fr(func(U)) for U in product([0, 1], repeat=size)}
        coeffs = {}
        for S in product([0, 1], repeat=size):
            tot = Fr(0)
            for U in product([0, 1], repeat=size):
                if all(u <= s for u, s in zip(U, S)):
                    tot += (-1) ** (sum(S) - sum(U)) * table[U]
            if tot: coeffs[S] = tot
        def phi(ins):
            out = ZPoly()
            for S, cf in coeffs.items():
                mono = ZPoly.lift(cf)
                for y, sbit in zip(ins, S):
                    if sbit: mono = mono * y
                out = out + mono
            return out
        return phi, (lambda xs: Poly.const(table[tuple(xs)]))
    scopes = []; facs = []; facs_b = []
    for v in A + B:
        inc = [e for e, (a, b) in enumerate(edges6) if v in (a, b)]
        sign = 1 if v in A else -1
        phi, fb = mle(len(inc), lambda U, sign=sign: int(sum(sign * (1 + u) for u in U) % 3 == 0))
        scopes.append(inc); facs.append(phi); facs_b.append(fb)
    results.append(run_case("nowhere-zero Z_3-flows of K_{2,3} over Q", n, scopes, facs, facs_b, 1, (3 + 3 + 2 + 2 + 2) * 2, 6))
    npass = sum(results)
    print(f"{npass} checks passed, {len(results) - npass} failed")
    sys.exit(0 if all(results) else 1)

if __name__ == "__main__":
    main()
