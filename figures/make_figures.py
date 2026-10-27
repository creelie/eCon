#!/usr/bin/env python3
"""Figures of the paper (PDF and PNG at 200 dpi). Run from any directory."""
import math, os, itertools, random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

HERE = os.path.dirname(os.path.abspath(__file__))

def save(fig, name):
    fig.savefig(os.path.join(HERE, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(HERE, name + ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------- figure 1
# The exponent of Theorem 6.7 as a function of the degree threshold b and of
# lambda = log n: leaves cost 2 ln(b+2)/(b+1) per vertex (natural log units),
# leaf cost (b+2) nu^2 / sqrt(lambda) per vertex (constant C_0 omitted).
def fig_branch():
    b = np.arange(2, 81)
    lam = 2.0 ** np.linspace(8, 24, 160)
    B, L = np.meshgrid(b, lam)
    nu = np.log2(L)
    leaves = 2 * np.log(B + 2) / (B + 1) / math.log(2)
    leaf = (B + 2) * nu ** 2 / np.sqrt(L)
    Z = leaves + leaf
    fig = plt.figure(figsize=(7.2, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(B, np.log2(L), np.log10(Z), cmap=cm.viridis, linewidth=0, antialiased=True, alpha=0.9)
    bl = np.floor(lam ** 0.25)
    nul = np.log2(lam)
    zl = 2 * np.log(bl + 2) / (bl + 1) / math.log(2) + (bl + 2) * nul ** 2 / np.sqrt(lam)
    ax.plot(bl, np.log2(lam), np.log10(zl) + 0.03, color="crimson", lw=2.5, zorder=10, label=r"$b=\lfloor\lambda^{1/4}\rfloor$")
    ax.set_xlabel(r"degree threshold $b$")
    ax.set_ylabel(r"$\log_2\lambda$")
    ax.set_zlabel(r"$\log_{10}$ of the exponent per vertex")
    ax.view_init(elev=32, azim=-140)
    ax.legend(loc="upper right")
    save(fig, "fig_branch")

# ---------------------------------------------------------------- figure 2
# The bisection profile q_{i,j} of Corollary 6.4 for a sparse graph on 14
# vertices, computed directly from the 2^14 subsets.
def fig_profile():
    rng = random.Random(2026)
    n = 14
    edges = set()
    for v in range(n):                       # a Hamilton cycle plus random chords
        edges.add((v, (v + 1) % n))
    while len(edges) < 21:
        u, v = rng.sample(range(n), 2)
        if (u, v) not in edges and (v, u) not in edges:
            edges.add((u, v))
    edges = sorted(edges)
    prof = np.zeros((n + 1, len(edges) + 1), dtype=np.int64)
    for mask in range(1 << n):
        i = bin(mask).count("1")
        j = sum(1 for u, v in edges if ((mask >> u) & 1) != ((mask >> v) & 1))
        prof[i, j] += 1
    fig = plt.figure(figsize=(7.4, 5.4))
    ax = fig.add_subplot(111, projection="3d")
    I, J = np.nonzero(prof)
    vals = prof[I, J]
    colors = cm.plasma(np.log1p(vals) / np.log1p(vals.max()))
    ax.bar3d(I - 0.4, J - 0.4, np.zeros_like(vals), 0.8, 0.8, np.log10(vals + 1), color=colors, shade=True)
    ax.set_xlabel(r"side size $i$")
    ax.set_ylabel(r"cut size $j$")
    ax.set_zlabel(r"$\log_{10}(1+q_{i,j})$")
    ax.view_init(elev=30, azim=-55)
    save(fig, "fig_profile")
    return n, len(edges), int(prof[:, :].sum()), int(np.max(np.nonzero(prof.sum(axis=0))[0])), int(np.max(np.nonzero(prof[n // 2])[0]))

# ---------------------------------------------------------------- figure 3
# The three costs of (3.x) with the parameters of Lemma 4.8, per variable,
# as a function of lambda and of the degree ratio kappa.
def fig_costs():
    beta = 2
    lam = 2.0 ** np.linspace(34, 60, 140)
    kap = np.arange(1, 21)
    K, L = np.meshgrid(kap, lam)
    nu = np.ceil(np.log2(L))
    r = np.floor(np.sqrt(L) / nu ** 2)
    Q = np.floor(L / (128 * nu))
    R = (beta + 3) * np.ceil(np.sqrt(L))
    tau = 1 / (r * L)
    logD = 2 * Q * np.log2(r)
    c1 = (K + 1) * (R / Q) * np.log2(2 * Q / R)
    c2 = tau * (logD + 1)
    c3 = 1 / r
    Z = c1 + c2 + c3
    Z = np.where((6 * R <= Q) & (r >= 2), Z, np.nan)
    fig = plt.figure(figsize=(7.2, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(K, np.log2(L), np.log10(Z), cmap=cm.cividis, linewidth=0, antialiased=True)
    ax.set_xlabel(r"degree ratio $\kappa$")
    ax.set_ylabel(r"$\log_2\lambda$")
    ax.set_zlabel(r"$\log_{10}$ of the exponent per variable")
    ax.view_init(elev=28, azim=-135)
    save(fig, "fig_costs")

# ---------------------------------------------------------------- figure 4
# The random-cluster partition function Z(G;q,v) of the Petersen graph over the
# (q,v) plane, from its 2^15 edge subsets: Z = sum_A q^{c(A)} v^{|A|}. The
# surface shows sign(Z) log10(1+|Z|); the lines are the hyperbolas
# (x-1)(y-1) = q at integer q >= 1 evaluated by Theorem 7.9 (fixed q, every v),
# the points the colourings (v=-1) and the nowhere-zero flows (v=-q).
def fig_tutte():
    n = 10
    outer = [(i, (i + 1) % 5) for i in range(5)]
    inner = [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    spokes = [(i, 5 + i) for i in range(5)]
    edges = outer + inner + spokes
    m = len(edges)
    # coefficients a[c][k] = number of edge subsets with k edges and c components
    coef = np.zeros((n + 1, m + 1), dtype=np.int64)
    for mask in range(1 << m):
        parent = list(range(n))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        k = 0
        for e in range(m):
            if (mask >> e) & 1:
                k += 1
                a, b = find(edges[e][0]), find(edges[e][1])
                if a != b:
                    parent[a] = b
        c = sum(1 for v in range(n) if find(v) == v)
        coef[c][k] += 1
    def Z(q, v):
        tot = 0.0
        for c in range(1, n + 1):
            for k in range(m + 1):
                if coef[c][k]:
                    tot += coef[c][k] * q ** c * v ** k
        return tot
    qs = np.linspace(0.5, 6.0, 111)
    vs = np.linspace(-6.0, 3.0, 181)
    Q, V = np.meshgrid(qs, vs)
    Zv = np.vectorize(Z)(Q, V)
    S = np.sign(Zv) * np.log10(1 + np.abs(Zv))
    fig = plt.figure(figsize=(7.4, 5.4))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(Q, V, S, cmap=cm.coolwarm, linewidth=0, antialiased=True, alpha=0.85)
    for q in range(1, 7):
        vv = np.linspace(-6.0, 3.0, 200)
        zz = np.array([Z(q, v) for v in vv])
        ax.plot(np.full_like(vv, q), vv, np.sign(zz) * np.log10(1 + np.abs(zz)) + 0.05, color="k", lw=1.6)
        for v, col, lab in ((-1.0, "crimson", "colourings"), (-float(q), "darkgreen", "flows")):
            z = Z(q, v)
            ax.scatter([q], [v], [np.sign(z) * np.log10(1 + abs(z)) + 0.1], color=col, s=28, zorder=20,
                       label=lab if q == 1 else None)
    ax.set_xlabel(r"$q$")
    ax.set_ylabel(r"$v$")
    ax.set_zlabel(r"$\mathrm{sign}(Z)\,\log_{10}(1+|Z|)$")
    ax.view_init(elev=30, azim=-125)
    ax.legend(loc="upper left")
    save(fig, "fig_tutte")
    return int(Z(3, -1)), int(Z(3, -3))

if __name__ == "__main__":
    fig_branch()
    print("profile graph:", fig_profile())
    fig_costs()
    print("Petersen P(3), F(3):", fig_tutte())
