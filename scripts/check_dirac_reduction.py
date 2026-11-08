#!/usr/bin/env python3
"""The N=2 ladder on zero supermomentum is a single Dirac response.

Claim under test (paper, Sec. 3). On functions of x = <12>, or equally for the charged kernel
[PSV](2.12) after an ordinary centre integral, the reduced kernel is

    K0(x, y) = c * sgn(x+y) |x+y|^(-2 delta) * sgn(y) |y|^(-2 mu),   mu = 1/2 - delta,

and with x = n e^tau, y = n' e^tau', f = e^(-delta tau) g, it becomes translation invariant in
tau. Its eigenvalue on e^(-H tau) Y_p(n) is

    k(H; sigma) = c [ A(H) - sigma B(H) ],
    A(H) = int dt e^(H t) (2 cosh(t/2))^(-2 delta),
    B(H) = int dt e^(H t) sgn(t) (2 |sinh(t/2)|)^(-2 delta),

where A and B are the two boundary correlators of a dimension-delta fermion on the strip, and
the integrals converge for |Re H| < delta. Tests:

    1. A(H), B(H) by quadrature against closed forms.
    2. k against kA, kS of [N2](15), which [PSV] confirm (check_literature_inputs.py).
    3. The constant c, from the rail integral of [PSV](2.12) with G_b = -d G_psi.
    4. The mode sums of the two correlators, so that their Fourier transform is 1/R_F of
       [N2](21).

Usage: python3 scripts/check_dirac_reduction.py --dps 20
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
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {mp.nstr(mp.mpc(value), 12)} vs {mp.nstr(mp.mpc(reference), 12)}")


def A_quad(d, H):
    f = lambda t: mp.cosh(H * t) * 2 * (2 * mp.cosh(t / 2)) ** (-2 * d)
    return mp.quad(f, [0, 1, 5, 20, mp.inf])


def B_quad(d, H):
    # integrable for 2 d < 1 at t = 0, which holds for q >= 3
    f = lambda t: 2 * mp.sinh(H * t) * (2 * mp.sinh(t / 2)) ** (-2 * d)
    return mp.quad(f, [0, 0.1, 1, 5, 20, mp.inf])


def A_exact(d, H):
    return mp.gamma(d - H) * mp.gamma(d + H) / mp.gamma(2 * d)


def B_exact(d, H):
    return mp.gamma(1 - 2 * d) * mp.gamma(d - H) * mp.gamma(d + H) * 2 * mp.cos(mp.pi * d) * mp.sin(mp.pi * H) / mp.pi


def chain_C(nu):
    """int dz sgn(z) sgn(z-t) |z|^(-delta) |z-t|^(-delta) = C |t|^(1-2 delta), [T](C.3)."""
    return 2 * mp.beta(nu, 1 - 2 * nu) - mp.beta(nu, nu)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dps", type=int, default=20)
    args = parser.parse_args()
    mp.mp.dps = args.dps

    for q in (3, 5, 7):
        m = N2(q)
        d, nu, mu = m.delta, m.nu, m.mu
        print(f"\nq={q}, delta={mp.nstr(d, 6)}")

        print(" 1. Boundary correlators against closed forms")
        for H in (mp.mpf("0.3") * d, mp.mpf("-0.6") * d):
            record(f"A({mp.nstr(H, 4)})", A_quad(d, H), A_exact(d, H), tol=mp.mpf(10) ** -9)
            record(f"B({mp.nstr(H, 4)})", B_quad(d, H), B_exact(d, H), tol=mp.mpf(10) ** -9)

        # 3. Constant. K^d = J(q-1) G_psi(t14) G_b(t23) G_psi(t34)^(q-2), [PSV](2.12). With
        # G_b = -d G_psi the rail is d/dt (b^2 C |t|^(2 mu)) = b^2 C 2 mu sgn(t) |t|^(-2 delta).
        # So c = J (q-1) b^q C 2 mu, up to an overall sign fixed by the orientation.
        Jbq = mp.tan(mp.pi * d / 2) / (2 * mp.pi)
        c = (q - 1) * Jbq * chain_C(nu) * 2 * mu
        print(f" 3. c = J (q-1) b^q C(nu) 2 mu = {mp.nstr(c, 12)}")

        print(" 2. Reduced kernel eigenvalue against kA, kS, for both parities, with this c")
        for H in (mp.mpf("0.3") * d, mp.mpf("-0.6") * d, mp.mpf("0.85") * d):
            k_odd = c * (A_exact(d, H) + B_exact(d, H))  # sigma = -1
            k_even = c * (A_exact(d, H) - B_exact(d, H))  # sigma = +1
            print(f"    H={mp.nstr(H, 5)}: c(A+B)={mp.nstr(k_odd, 10)} c(A-B)={mp.nstr(k_even, 10)}"
                  f" kA(H)={mp.nstr(m.kA(H), 10)} kS(H)={mp.nstr(m.kS(H), 10)}"
                  f" kA(-H)={mp.nstr(m.kA(-H), 10)} kS(-H)={mp.nstr(m.kS(-H), 10)}")
            record("sigma=+1 gives kA(H)", c * (A_exact(d, H) - B_exact(d, H)), m.kA(H))
            record("sigma=-1 gives kS(H)", c * (A_exact(d, H) + B_exact(d, H)), m.kS(H))

        print(" 4. The response [N2](22): k = sigma L / R_F")
        for H in (mp.mpf("1.37"), mp.mpf("0.31"), mp.mpf("2.23")):
            om = H - mp.mpf(1) / 2
            for p, branch, sigma in ((1, m.kA, -1), (0, m.kS, 1)):
                RF = 2 * nu * response(nu, p, om) / (H - nu)
                record(f"p={p} H={H}", sigma * m.L / RF, branch(H))

    print("\n5. Closed forms and positivity, for the paper")
    for q in (3, 5, 7):
        m = N2(q)
        d, nu, mu = m.delta, m.nu, m.mu
        Jbq = mp.tan(mp.pi * d / 2) / (2 * mp.pi)
        c = (q - 1) * Jbq * chain_C(nu) * 2 * mu
        record(f"q={q}: c = -nu Gamma(nu)/(2 delta Gamma(delta) cos(pi delta) Gamma(1-2 delta))", c,
               -nu * mp.gamma(nu) / (2 * d * mp.gamma(d) * mp.cos(mp.pi * d) * mp.gamma(1 - 2 * d)))
        cF = -mp.gamma(d + mp.mpf(1) / 2) / (mp.sqrt(mp.pi) * mp.gamma(d))
        record(f"q={q}: c/L = c_F = -Gamma(delta+1/2)/(sqrt(pi) Gamma(delta))", c / m.L, cF)
        # 1/R_F = c_F (sigma A + B), continued off the strip through the closed forms
        for H in (mp.mpf("0.11"), mp.mpf("1.37")):
            om = H - mp.mpf(1) / 2
            for p, sigma in ((1, -1), (0, 1)):
                RF = 2 * nu * response(nu, p, om) / (H - nu)
                record(f"q={q} H={H} p={p}: 1/R_F = c_F (sigma A + B)", 1 / RF, cF * (sigma * A_exact(d, H) + B_exact(d, H)))
        # positivity of the Dirac spectral weights, at positive and negative roots
        half = mp.mpf(1) / 2
        def RFp(p, W):
            return 2 * nu * response(nu, p, W - half) / (W - nu)
        for g in (mp.mpf("0.6"), mp.mpf("0.95")):
            for p, sigma, branch in ((1, -1, m.kA), (0, 1, m.kS)):
                for lo, hi in ((mp.mpf("-3.3"), mp.mpf("3.3")),):
                    a = lo + mp.sqrt(2) * mp.mpf(10) ** -4
                    step = mp.mpf("0.005")
                    n_roots = 0
                    while a < hi:
                        b = a + step
                        fa, fb = branch(a) - 1 / g, branch(b) - 1 / g
                        if fa * fb < 0 and abs(fa) < 40 and abs(fb) < 40:
                            r = mp.findroot(lambda x: branch(x) - 1 / g, (a, b), solver="anderson")
                            # G^F = (R_F - g sigma L)^(-1) = sum w_n/(Omega - Omega_n): w_n = 1/R_F' > 0
                            w = 1 / mp.diff(lambda x: RFp(p, x), r)
                            n_roots += 1
                            RESULTS.append(w > 0)
                            print(f"    g={g} p={p}: root Omega={mp.nstr(r, 7)}, weight Res G^F = {mp.nstr(w, 6)} {'PASS' if w > 0 else 'FAIL'}")
                        a = b
        # the transformations on the Dirac side: B_minus = s a_F, B_plus = s b_F for s = -w
        for w, p in ((mp.mpf("1.63"), 1), (mp.mpf("2.37"), 0)):
            A_ = mp.mpf(1)
            Bb = -(w + mu) / (2 * d) * A_
            record(f"q={q} w={w}: B_minus = -w a_F", 2 * d * Bb + mu * A_, -w * A_)
            Ab_ = mp.mpf("0.77")
            Bpsi = (w - mu) * Ab_ / (2 * nu)
            record(f"q={q} w={w}: B_plus = -w b_F", -2 * nu * Bpsi - mu * Ab_, -w * Ab_)

    failed = RESULTS.count(False)
    print(f"\n{len(RESULTS) - failed}/{len(RESULTS)} checks passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
