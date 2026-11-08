#!/usr/bin/env python3
"""Off-shell Ward identity, spectrum multiplets, BPS threshold and vacuum energy; data for figures.

Claims under test (paper, Sec. 5):

  1. The neutral 2x2 resolvent is two Dirac resolvents, at every complex w and real g:
       G_{g,p}(w) M_p = P_-(w) sigma_p L G^F_{g,p}(w + 1/2)
                        + P_+(w) sigma_{1-p} L G^F_{g,1-p}(w - 1/2),
     with P_pm the projectors on u_pm(w) = (1, (+-w - mu)/(2 delta)) of [N2](14).
     The second term uses the second form of R_F, [N2] below (23).
  2. Quadratic form: X^T M_p X = 2 sigma L A (mu A + 2 delta B), the boundary term being A times
     Q X with Q X = (d_tau + mu) A + 2 delta B, since A d_tau A is a total derivative.
  3. Every root organizes into energies (H-1/2, H, H, H+1/2) for the two scalars and the Dirac
     field. The su(1|1) central charge 2 omega = 2H - 1 vanishes at H = 1/2, which is g = g_*.
  4. Vacuum energy per multiplet: (|H-1/2| + (H+1/2))/2 - H, zero for H >= 1/2.

With --write it records H_S(g), H_A(g) for q=3,5 in results/data/paper/ for the figure.

Usage: python3 scripts/check_ward_and_vacuum.py --dps 20 [--write]
"""

import argparse
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_robin_identities import N2, response  # noqa: E402

RESULTS = []


def record(label, value, reference, tol=None):
    tol = tol if tol is not None else mp.mpf(10) ** (-mp.mp.dps + 8)
    diff = abs(mp.mpc(value) - mp.mpc(reference))
    ok = diff / max(1, abs(mp.mpc(reference))) < tol
    RESULTS.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {mp.nstr(mp.mpc(value), 10)} vs {mp.nstr(mp.mpc(reference), 10)}")


def RF_first(m, p, Omega):
    """[N2](21)."""
    return 2 * m.nu * response(m.nu, p, Omega - mp.mpf(1) / 2) / (Omega - m.nu)


def RF_second(m, p, Omega):
    """[N2], the form below (23)."""
    return 2 * m.nu * response(m.nu, 1 - p, Omega + mp.mpf(1) / 2) / (Omega + mp.mpf(1) / 2 + m.mu)


def check_ward(m, g):
    half = mp.mpf(1) / 2
    for w in (mp.mpf("0.83"), mp.mpc("1.3", "0.7"), mp.mpc(0, "1.9")):
        for p in (0, 1):
            sig = m.sigma(p)
            D, M = m.Dmat(p, w), m.Mmat(p)
            lhs = mp.inverse(mp.inverse(D) - g * M) * M
            up = mp.matrix([1, (w - m.mu) / (2 * m.delta)])
            um = mp.matrix([1, (-w - m.mu) / (2 * m.delta)])
            U = mp.matrix([[um[0], up[0]], [um[1], up[1]]])
            Ui = mp.inverse(U)
            Pm = um * mp.matrix([[Ui[0, 0], Ui[0, 1]]])
            Pp = up * mp.matrix([[Ui[1, 0], Ui[1, 1]]])
            Gm = 1 / (RF_first(m, p, w + half) - g * sig * m.L)
            Gp = 1 / (RF_second(m, 1 - p, w - half) - g * (-sig) * m.L)
            rhs = Pm * (sig * m.L * Gm) + Pp * ((-sig) * m.L * Gp)
            for i in range(2):
                for j in range(2):
                    record(f"g={g} w={mp.nstr(w, 3)} p={p} entry {i}{j}", lhs[i, j], rhs[i, j])


def roots(m, k, g, lo, hi, n=800):
    out = []
    # an irrational offset keeps the grid off the exact Gamma-function poles at H = delta + n
    grid = [lo + mp.sqrt(2) * mp.mpf(10) ** -5 + (hi - lo) * mp.mpf(i) / n for i in range(n + 1)]
    for a, b in zip(grid, grid[1:]):
        fa, fb = k(a) - 1 / g, k(b) - 1 / g
        if fa * fb < 0 and abs(fa) < 30 and abs(fb) < 30:
            out.append(mp.findroot(lambda x: k(x) - 1 / g, (a, b), solver="anderson"))
    return out


def check_threshold(m):
    half = mp.mpf(1) / 2
    gstar = 1 / m.kS(half)
    print(f"  g_* = 1/kS(1/2) = {mp.nstr(gstar, 10)}")
    below = roots(m, m.kS, gstar * mp.mpf("0.6"), m.delta + mp.mpf("0.002"), mp.mpf("0.999"))
    above = roots(m, m.kS, min(gstar * mp.mpf("1.6"), mp.mpf(1)), m.delta + mp.mpf("0.002"), mp.mpf("0.999"))
    record("below g_*, lowest H_S is under 1/2", 1 if below[0] < half else 0, 1)
    record("above g_*, lowest H_S is over 1/2", 1 if above[0] > half else 0, 1)
    for label, rts in (("below", below), ("above", above)):
        H = rts[0]
        E = (abs(H - half) + H + half) / 2 - H
        print(f"    {label} g_*: H_S = {mp.nstr(H, 8)}, vacuum energy of the multiplet = {mp.nstr(E, 8)}")
    H = below[0]
    record("below g_*: E = 1/2 - H", (abs(H - half) + H + half) / 2 - H, half - H)
    H = above[0]
    record("above g_*: E = 0", (abs(H - half) + H + half) / 2 - H, 0)
    # the charged spectrum is exactly the neutral one shifted: test against the 2x2 determinant
    for g in (mp.mpf("0.7"), mp.mpf("0.95")):
        for Hr in roots(m, m.kA, g, mp.mpf("1.2"), mp.mpf(4)):
            for w, p in ((Hr - half, 1), (Hr + half, 0)):
                kb = m.kbar(p, w)
                det = (1 - g * kb[0, 0]) * (1 - g * kb[1, 1]) - g**2 * kb[0, 1] * kb[1, 0]
                record(f"A root {mp.nstr(Hr, 6)}: det at w={mp.nstr(w, 6)}, p={p}", det, 0, tol=mp.mpf(10) ** -8)


def write_data(root):
    os.makedirs(root, exist_ok=True)
    for q in (3, 5):
        m = N2(q)
        rows = []
        for i in range(1, 121):
            g = mp.mpf(i) / 120
            s = roots(m, m.kS, g, m.delta + mp.mpf("0.0005"), mp.mpf("1.0005"), n=600)
            a = roots(m, m.kA, g, mp.mpf("1.0005"), mp.mpf("1.6") + m.delta, n=600)
            s2 = roots(m, m.kS, g, mp.mpf("1.5"), mp.mpf("3.2"), n=600)
            if s and a and s2:
                rows.append((g, s[0], a[0], s2[0]))
        path = os.path.join(root, f"spectrum_q{q}.dat")
        with open(path, "w") as f:
            f.write("# g H_S0 H_A0 H_S1 ; lowest roots of k_alpha(H)=1/g; scripts/check_ward_and_vacuum.py\n")
            f.write("g HS0 HA0 HS1\n")
            for r in rows:
                f.write(" ".join(mp.nstr(x, 8) for x in r) + "\n")
        print(f"  wrote {path}, {len(rows)} rows")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dps", type=int, default=20)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    mp.mp.dps = args.dps
    for q in (3, 5, 7):
        m = N2(q)
        print(f"\nq={q}")
        for g in (mp.mpf("0.6"), mp.mpf(1)):
            check_ward(m, g)
        # the quadratic form
        A, B = mp.mpf("0.7"), mp.mpf("-1.3")
        for p in (0, 1):
            M = m.Mmat(p)
            X = mp.matrix([A, B])
            record(f"X^T M X, p={p}", (X.T * M * X)[0], 2 * m.sigma(p) * m.L * A * (m.mu * A + 2 * m.delta * B))
        check_threshold(m)
    if args.write:
        write_data(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "data", "paper"))
    failed = RESULTS.count(False)
    print(f"\n{len(RESULTS) - failed}/{len(RESULTS)} checks passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
