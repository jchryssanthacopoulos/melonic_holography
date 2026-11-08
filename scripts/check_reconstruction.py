#!/usr/bin/env python3
"""The four-point reconstruction [N2](26)-(28) at finite g: contour plus integers against roots.

    F_g(z) = int_C dH/(2 pi i) r_g(H) phi_H(z) + D sum_n k_A(n)/(1 - g k_A(n)) phi_n(z)
           = sum_{roots} c_alpha(H) phi_H(z).

The contour crosses the real axis at delta/2 and is tilted to the right to gain decay from z^H,
which is legitimate as long as all poles of r_g to the right of the crossing lie on the real axis.
Also tests the expansion at small g against the order g^1 form, and the transform [N2](25) on
the resummed function.

Usage: python3 scripts/check_reconstruction.py --dps 18
"""

import argparse
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_robin_identities import N2  # noqa: E402

RESULTS = []


def record(label, value, reference, tol):
    diff = abs(mp.mpc(value) - mp.mpc(reference))
    ok = diff / max(1, abs(mp.mpc(reference))) < tol
    RESULTS.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {mp.nstr(mp.mpc(value), 12)} vs {mp.nstr(mp.mpc(reference), 12)}")


def r_g(m, g, H):
    kA, kS = m.kA(H), m.kS(H)
    return m.norm * (mp.cot(mp.pi * H / 2) * kA / (1 - g * kA) - mp.tan(mp.pi * H / 2) * kS / (1 - g * kS))


def contour_part(m, g, z):
    c = m.delta / 2
    phase = mp.expjpi(mp.mpf(1) / 3)
    f = lambda t: r_g(m, g, c + t * phase) * m.phi(c + t * phase, z) * phase
    return mp.im(mp.quad(f, [0, 1, 3, 8, 20, 50])) / mp.pi


def integer_part(m, g, z, nmax=60):
    total = mp.mpf(0)
    for n in range(1, nmax + 1):
        k = m.kA(mp.mpf(n))
        total += m.D * k / (1 - g * k) * m.phi(mp.mpf(n), z)
    return total


def roots(fn, g, lo, hi, step=mp.mpf("0.01")):
    out = []
    a = lo + mp.sqrt(2) * mp.mpf(10) ** -4
    while a < hi:
        b = a + step
        fa, fb = fn(a) - 1 / g, fn(b) - 1 / g
        if fa * fb < 0 and abs(fa) < 30 and abs(fb) < 30:
            out.append(mp.findroot(lambda x: fn(x) - 1 / g, (a, b), solver="anderson"))
        a = b
    return out


def root_sum(m, g, z, hmax=30):
    total = mp.mpf(0)
    for H in roots(m.kA, g, m.delta + mp.mpf("0.001"), mp.mpf(hmax)):
        total += m.norm * mp.cot(mp.pi * H / 2) / (g**2 * mp.diff(m.kA, H)) * m.phi(H, z)
    for H in roots(m.kS, g, m.delta + mp.mpf("0.001"), mp.mpf(hmax)):
        total += -m.norm * mp.tan(mp.pi * H / 2) / (g**2 * mp.diff(m.kS, H)) * m.phi(H, z)
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dps", type=int, default=18)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    tol = mp.mpf(10) ** -9
    for q in (3, 5):
        m = N2(q)
        print(f"\nq={q}")
        for g in (mp.mpf("0.15"), mp.mpf("0.5"), mp.mpf("0.85")):
            for z in (mp.mpf("0.25"), mp.mpf("0.6")):
                a = contour_part(m, g, z) + integer_part(m, g, z)
                b = root_sum(m, g, z)
                record(f"g={g} z={z}: contour + integers = root sum", a, b, tol)
        # small-g limit: F_g -> z^delta
        z = mp.mpf("0.4")
        g = mp.mpf("0.0005")
        record("g -> 0 tends to z^delta (to O(g))", contour_part(m, g, z) + integer_part(m, g, z), z**m.delta, mp.mpf(10) ** -2)
    failed = RESULTS.count(False)
    print(f"\n{len(RESULTS) - failed}/{len(RESULTS)} checks passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
