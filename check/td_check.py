"""Recompute H/I and min_i E_i/Inf_i for the affine transversal family TD(w,Q), Q prime.

Variables (x, y) with x in {0..w-1}, y in F_Q; terms t_{a,b} = {(x, a x + b) : x in {0..w-1}} for a, b in F_Q.
Exact integer Walsh transform of f = 1 - 2F, then logarithms. Usage: python3 td_check.py w Q
"""
import json
import sys

import numpy as np


def wht(v):
    v = v.astype(np.int64).copy()
    n = v.size
    h = 1
    while h < n:
        v = v.reshape(-1, 2, h)
        a = v[:, 0, :].copy()
        b = v[:, 1, :]
        v = np.stack([a + b, a - b], axis=1).reshape(n)
        h *= 2
    return v


def phi(p, q):
    s = p + q
    out = np.zeros_like(s)
    m = p > 0
    out[m] += p[m] * np.log2(s[m] / p[m])
    m = q > 0
    out[m] += q[m] * np.log2(s[m] / q[m])
    return out


def main():
    w, Q = int(sys.argv[1]), int(sys.argv[2])
    n = w * Q
    var = {(x, y): x * Q + y for x in range(w) for y in range(Q)}
    terms = []
    for a in range(Q):
        for b in range(Q):
            terms.append(sum(1 << var[(x, (a * x + b) % Q)] for x in range(w)))
    for s in range(len(terms)):
        for t in range(s + 1, len(terms)):
            assert bin(terms[s] & terms[t]).count("1") <= 1
    N = 1 << n
    idx = np.arange(N, dtype=np.int64)
    F = np.zeros(N, dtype=bool)
    for t in terms:
        F |= (idx & t) == t
    W = wht(1 - 2 * F.astype(np.int64))
    P = (W.astype(np.float64) ** 2) / float(N) ** 2
    nz = P > 0
    H = float(-(P[nz] * np.log2(P[nz])).sum())
    infs, ratios = [], []
    for i in range(n):
        lo = idx[(idx >> i) & 1 == 0]
        p, q = P[lo], P[lo | (1 << i)]
        inf = float(q.sum())
        infs.append(inf)
        ratios.append(float(phi(p, q).sum()) / inf)
    I = sum(infs)
    out = {"w": w, "Q": Q, "n": n, "terms": len(terms), "P_F0": float(1 - F.mean()), "H": H, "I": I,
           "H_over_I": H / I, "min_E_over_Inf": min(ratios), "max_E_over_Inf": max(ratios)}
    print(json.dumps(out))


if __name__ == "__main__":
    main()
