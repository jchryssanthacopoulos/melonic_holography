#!/usr/bin/env python3
"""The super-centre integral of the N=2 chiral-antichiral four-point function.

Thong reduces a ladder by integrating one pair's centre Y. Here pair (12) is integrated over its
whole supertranslation orbit (Y, eta, etabar), with pair (34) fixed, and the result is compared
with the Dirac propagator of the paper. Conventions: research/tex/n2_dirac_ladders/, with
z = <12><34>/(<14><32>) (the chi of Bulycheva) and chi = z/(z-1) (the chi of [PSV]); the
normalized correlator is Fcal(z) = -F_PSV(chi); the Berezin measure has int d^2eta eta etabar = 1.

  1. Grassmann algebra: the generators G_{-1/2}, Gbar_{-1/2} move the pair along an orbit on which
     <12> = x, <34> = y, <14> = Y + s + eta etabar, <32> = s - Y + eta etabar, s = (x+y)/2. So the
     eta etabar part of z is -2 s z0^2/(xy), with z0 the body.
  2. Counting: on four chiral-antichiral points the algebra acts with rank 3|4, so z is the only
     invariant and there is no nilpotent one.
  3. The seed Fcal_0 = sgn(z)|z|^a: direct quadrature against Beta functions for -1/2 < a < 0, where
     the integral converges; continued to a = delta, it is 8 pi N_delta L_q times the Dirac
     correlator with its second argument reflected.
  4. One superfamily, assembled from the component families of [T](C.4) and the even-sector
     analogue: the super-centre integral is 2 d_t of the ordinary one divided by sinh(t/2) or
     cosh(t/2), and it collapses to one exponential e^{-H|t|}, weight -4 pi tan(pi H/2) (A type) or
     4 pi cot(pi H/2) (S type). Also for H < 1/2, where the ordinary integral diverges.
  5. Finite g from [PSV](4.10), (4.11), (4.26): at each root, the super-centre weight of the
     literature's four-point function, -4 Res[(2h-1) k0^{psi,s}/D^a] or 4 Res[(2h-1) k0^{psi,a}/D^s],
     equals the Dirac weight -/+ 4 pi N_delta L_q w_n, w_n = 1/R_F'(H_n). Includes the exceptional
     roots at g = 1 and roots below g_*.
  6. The mechanism, from the Ward identity: G_{g,p,11}(w) = -[G^F_{g,p}(w+1/2) + G^F_{g,1-p}(w-1/2)]/2w.

Usage: python3 scripts/check_super_centre.py --dps 20 [--parts 123456]
Exits non-zero if a check fails. See docs/derivations.md Sec. 12.
"""

import argparse
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from check_literature_inputs import PSV, _block_a, _block_s, _res  # noqa: E402
from check_robin_identities import N2  # noqa: E402

RESULTS = []
half = mp.mpf(1) / 2


def record(label, value, reference, tol=None):
    tol = tol if tol is not None else mp.mpf(10) ** (-mp.mp.dps + 8)
    diff = abs(mp.mpc(value) - mp.mpc(reference))
    ok = diff / max(1, abs(mp.mpc(reference))) < tol
    RESULTS.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {mp.nstr(mp.mpc(value), 10)} vs {mp.nstr(mp.mpc(reference), 10)}")
    return ok


def info(label, value):
    print(f"  [INFO] {label}: {value if isinstance(value, str) else mp.nstr(value, 10)}")


def Nn(h):
    return mp.gamma(h) ** 2 / (2 * mp.gamma(2 * h))


# ------------------------------------------------------------------ 1
def part_1():
    print("\n1. Brackets along the supertranslation orbit of the pair (12)")
    import sympy as sp

    from check_supercasimir import generators, lift, odd

    Y, x, y, t1, t2 = sp.symbols("Y x y tau1 tau2", real=True)
    eta, etab = odd("t3"), odd("tb3")

    def bracket(a, b):
        # <ij> = tau_ij - 2 thetabar_i theta_j - theta_i thetabar_i - theta_j thetabar_j, chiral i, antichiral j
        (ta, tha, thba), (tb, thb, thbb) = a, b
        return lift(ta) - lift(tb) - 2 * (thba * thb) - tha * thba - thb * thbb

    # The flow of X = eta G_{-1/2} + etabar Gbar_{-1/2}, from the generators of the paper
    gens = generators(1, sp.Integer(0), 0, 1)
    X = lambda f: eta * gens["Gm"](f) + etab * gens["Gbm"](f)
    th1, thb1 = odd("t1"), odd("tb1")
    tau = lift(sp.Symbol("tau1", real=True))
    ok = (X(tau) - (-(eta * th1) - etab * thb1)).is_zero() and (X(th1) - etab).is_zero() and (X(thb1) - eta).is_zero()
    ok = ok and X(X(tau)).is_zero()
    RESULTS.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] X tau = -eta theta - etabar thetabar, X theta = etabar, X thetabar = eta, X^2 tau = 0")

    def flow(t, th, thb):
        return lift(t) - eta * th - etab * thb, th + etab, thb + eta

    # invariance of a chiral-antichiral bracket under the flow of both points
    p1 = (t1, odd("t1"), odd("tb1"))
    p2 = (t2, odd("t2"), odd("tb2"))
    ok = (bracket(flow(*p1), flow(*p2)) - bracket(p1, p2)).is_zero()
    RESULTS.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] <12> is invariant under the supertranslation of both points")

    zero = lift(0)
    P1 = flow(Y + x / 2, zero, zero)
    P2 = flow(Y - x / 2, zero, zero)
    P3 = (y / 2, zero, zero)
    P4 = (-y / 2, zero, zero)
    s = (x + y) / 2
    checks = {
        "<12> = x": (bracket(P1, P2), lift(x)),
        "<34> = y": (bracket(P3, P4), lift(y)),
        "<14> = Y + s + eta etabar": (bracket(P1, P4), lift(Y + s) + eta * etab),
        "<32> = s - Y + eta etabar": (bracket(P3, P2), lift(s - Y) + eta * etab),
    }
    for label, (a, b) in checks.items():
        ok = (a - b).expand().is_zero()
        RESULTS.append(bool(ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}")

    def inverse(g):
        body = g.terms.get((), 0)
        nil = g - lift(body)
        out, term = lift(1 / body), lift(1)
        for _ in range(3):
            term = term * nil * lift(-1 / body)
            out = out + term * lift(1 / body)
        return out

    z = bracket(P1, P2) * bracket(P3, P4) * inverse(bracket(P1, P4)) * inverse(bracket(P3, P2))
    z0 = x * y / (s**2 - Y**2)
    top = z.terms.get(("t3", "tb3"), 0)
    ok = sp.simplify(z.terms.get((), 0) - z0) == 0 and sp.simplify(top - (-2 * s * z0**2 / (x * y))) == 0
    RESULTS.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] z = z0 - 2 s z0^2 eta etabar/(xy), z0 = xy/(s^2 - Y^2)")


# ------------------------------------------------------------------ 2
def part_2():
    print("\n2. Counting: rank of the algebra on four chiral-antichiral points, at a generic body")
    T = [mp.mpf("0.31"), mp.mpf("-1.7"), mp.mpf("2.2"), mp.mpf("0.9")]
    # even generators L_1, L_0, L_{-1} on (tau_1..tau_4); J_0 moves only odd coordinates
    even = mp.matrix([[1, 1, 1, 1], T, [t**2 for t in T]])
    # odd coordinates: thetabar of the chiral points 1, 3, theta of the antichiral points 2, 4.
    # G_{-1/2}, G_{+1/2} shift thetabar by 1 and tau; Gbar_{-1/2}, Gbar_{+1/2} shift theta.
    oddm = mp.matrix([[1, 1, 0, 0], [T[0], T[2], 0, 0], [0, 0, 1, 1], [0, 0, T[1], T[3]]])
    import sympy as sp

    re = sp.Matrix(3, 4, lambda i, j: sp.nsimplify(str(even[i, j])))
    ro = sp.Matrix(4, 4, lambda i, j: sp.nsimplify(str(oddm[i, j])))
    record("rank of the even block (3 means one even invariant)", re.rank(), 3)
    record("rank of the odd block (4 means no odd invariant)", ro.rank(), 4)


# ------------------------------------------------------------------ 3
def S_beta(a):
    return mp.beta(half, -a) + mp.beta(half + a, -a)


def seed_direct(a, C):
    """int_R du a |z|^(1+a), z = +-1/(C^2 - u^2), by u = C sin(phi) and u = C cosh(psi)."""
    with mp.workdps(mp.mp.dps + 30):
        inner = mp.quad(lambda ph: mp.cos(ph) ** (-1 - 2 * a), [0, mp.pi / 4, mp.pi / 2])
        outer = mp.quad(lambda ps: mp.sinh(ps) ** (-1 - 2 * a), [0, 1, mp.inf])
        return 2 * a * C ** (-1 - 2 * a) * (inner + outer)


def part_3():
    print("\n3. The seed sgn(z)|z|^a: direct integral, Beta functions, and the Dirac correlator")
    for a in (mp.mpf("-0.2"), mp.mpf("-0.35")):
        for C in (mp.mpf("1.3"), mp.mpf("0.45")):
            record(f"a={a} C={C}: int du a|z|^(1+a) = a C^(-1-2a) [B(1/2,-a) + B(1/2+a,-a)]", seed_direct(a, C), a * C ** (-1 - 2 * a) * S_beta(a), tol=mp.mpf(10) ** -12)
    for q in (3, 5, 7):
        m = N2(q)
        d = m.delta
        cfrak = mp.gamma(d + half) / (mp.sqrt(mp.pi) * mp.gamma(d))
        K = 8 * mp.pi * m.norm * m.L * cfrak
        for t in (mp.mpf("0.6"), mp.mpf("-1.1"), mp.mpf("2.3")):
            c, sh = mp.cosh(t / 2), mp.sinh(t / 2)
            # R/n: same sign = -2 cosh(t/2) J_c, opposite sign = -2 sinh(t/2) J_s, continued to a = delta
            same = -2 * d * S_beta(d) * c ** (-2 * d)
            opp = -2 * d * S_beta(d) * mp.sign(t) * abs(sh) ** (-2 * d)
            record(f"q={q} t={t}: n'=n, against 8 pi N L cfrak (2cosh(t/2))^(-2delta)", same, K * (2 * c) ** (-2 * d))
            record(f"q={q} t={t}: n'=-n, against 8 pi N L cfrak sgn t |2sinh(t/2)|^(-2delta)", opp, K * mp.sign(t) * abs(2 * sh) ** (-2 * d))
        b = mp.power(mp.tan(mp.pi * d / 2) / (2 * mp.pi), d)
        ell = 2 * mp.power(4, -m.nu) / b ** (q - 2)
        record(f"q={q}: 8 pi N_delta L_q b^2 = 2 l_psi in [PSV]'s normalization", 8 * mp.pi * m.norm * m.L * b**2, 2 * ell)


# ------------------------------------------------------------------ 4
def family(kind, H):
    if kind == "A":  # bottom in p = 1, top in p = 0
        return lambda chi: 2 * Nn(H) * _block_s(H, chi) + 2 * Nn(H + 1) * _block_a(H + 1, chi)
    return lambda chi: -2 * Nn(H) * _block_a(H, chi) - 2 * Nn(H + 1) * _block_s(H + 1, chi)


def P_ord(fun, t, same, U0=None):
    """Ordinary centre integral, chi = 1/(u^2 - sinh^2) (same sign) or 1/(cosh^2 - u^2) (opposite).

    With U0, integrate only |u| < U0; otherwise over all u."""
    s, c = abs(mp.sinh(t / 2)), mp.cosh(t / 2)
    f = (lambda u: fun(1 / (u**2 - s**2))) if same else (lambda u: fun(1 / (c**2 - u**2)))
    total = mp.quad(f, [0, s]) + mp.quad(f, [s, c])
    if U0 is None:
        total += mp.quad(f, [c, 2 * c, 10 * c, mp.inf])
    else:
        total += mp.quad(f, [c, (c + U0) / 2, U0])
    return 2 * total


def d_dt(F, t, e=mp.mpf("0.002")):
    return (F(t - 2 * e) - 8 * F(t - e) + 8 * F(t + e) - F(t + 2 * e)) / (12 * e)


def super_from_ordinary(fun, t, same, U0=None, tail=None):
    """R/n from the ordinary centre integral: 2 d_t P / sinh(t/2) or -2 d_t P / cosh(t/2)."""
    dP = d_dt(lambda tt: P_ord(fun, tt, same, U0), t)
    if tail is not None:
        dP += tail
    return 2 * dP / mp.sinh(t / 2) if same else -2 * dP / mp.cosh(t / 2)


def prediction(kind, H, t, same):
    if kind == "A":
        w = -4 * mp.pi * mp.tan(mp.pi * H / 2)
        return w * mp.exp(-H * abs(t)) * (1 if same else -mp.sign(t))
    w = 4 * mp.pi * mp.cot(mp.pi * H / 2)
    return w * mp.exp(-H * abs(t)) * (1 if same else mp.sign(t))


def part_4():
    print("\n4. One superfamily: the super-centre integral is a single exponential e^{-H|t|}")
    print(" R/n with n' = n and n' = -n; predicted -4 pi tan(pi H/2) [theta(t) n' + theta(-t) n] e^{-H|t|} (A),")
    print(" 4 pi cot(pi H/2) [theta(t) n + theta(-t) n'] e^{-H|t|} (S)")
    for kind in ("A", "S"):
        for H in (mp.mpf("1.3"), mp.mpf("1.7"), mp.mpf("2.4")):
            fun = family(kind, H)
            for t in (mp.mpf("0.7"), mp.mpf("-0.9")):
                for same in (True, False):
                    val = super_from_ordinary(fun, t, same)
                    record(f"{kind} H={H} t={t} {'same' if same else 'opposite'} sign", val, prediction(kind, H, t, same), tol=mp.mpf(10) ** -8)
            # the ordinary centre integral of the same superfamily has two exponentials
            Pa, Pb = P_ord(fun, mp.mpf("0.7"), True), P_ord(fun, mp.mpf("1.4"), True)
            info(f"   ordinary P at t=0.7, 1.4 (two exponentials e^(-(H-1/2)|t|), e^(-(H+1/2)|t|))", f"{mp.nstr(Pa, 8)}, {mp.nstr(Pb, 8)}")
    print(" Below H = 1/2, where the ordinary centre integral of the bottom component diverges")
    U0 = mp.mpf(5)
    for kind, H in (("S", mp.mpf("0.3")), ("A", mp.mpf("0.4")), ("S", mp.mpf("0.45"))):
        fun = family(kind, H)
        dfun = lambda chi: mp.diff(fun, chi)
        for t in (mp.mpf("0.7"), mp.mpf("-0.9")):
            for same in (True, False):
                s, c = abs(mp.sinh(t / 2)), mp.cosh(t / 2)
                if same:
                    # d_t chi = chi^2 sinh(t)/2 for chi = 1/(u^2 - s^2), which lies in (0,1) for u > U0
                    g = lambda u: dfun(1 / (u**2 - s**2)) * (1 / (u**2 - s**2)) ** 2 * mp.sinh(t) / 2
                else:
                    g = lambda u: -dfun(1 / (c**2 - u**2)) * (1 / (c**2 - u**2)) ** 2 * mp.sinh(t) / 2
                tail = 2 * mp.quad(g, [U0, 2 * U0, 10 * U0, mp.inf])
                val = super_from_ordinary(fun, t, same, U0=U0, tail=tail)
                record(f"{kind} H={H} t={t} {'same' if same else 'opposite'} sign", val, prediction(kind, H, t, same), tol=mp.mpf(10) ** -7)
        grow = [P_ord(fun, mp.mpf("0.7"), True, U0=U) for U in (mp.mpf(50), mp.mpf(500))]
        info(f"   {kind} H={H}: ordinary P cut at |u| = 50, 500 (it diverges)", f"{mp.nstr(grow[0], 8)}, {mp.nstr(grow[1], 8)}")


# ------------------------------------------------------------------ 5
def RF(m, p, W):
    return 2 * m.nu * m.Rnu(p, W - half) / (W - m.nu)


def scan_roots(k, g, lo, hi):
    out = []
    x = lo + mp.sqrt(2) * mp.mpf(10) ** -4
    step = mp.mpf("0.01")
    while x < hi:
        fa, fb = k(x) - 1 / g, k(x + step) - 1 / g
        if fa * fb < 0 and abs(fa) < 30 and abs(fb) < 30:
            out.append(mp.findroot(lambda yy: k(yy) - 1 / g, (x, x + step), solver="anderson"))
        x += step
    return out


def part_5():
    print("\n5. Finite g: the super-centre weights of [PSV]'s four-point function against the Dirac weights")
    for q, gs in ((3, ("0.15", "0.5", "0.9", "1")), (5, ("0.3", "0.7", "1"))):
        m = N2(q)
        psv = PSV(q)
        record(f"q={q}: k_A(H) = -L_q / R_F,1(H) at H=0.83", m.kA(mp.mpf("0.83")), -m.L / RF(m, 1, mp.mpf("0.83")), tol=mp.mpf(10) ** -12)
        record(f"q={q}: k_S(H) = +L_q / R_F,0(H) at H=0.83", m.kS(mp.mpf("0.83")), m.L / RF(m, 0, mp.mpf("0.83")), tol=mp.mpf(10) ** -12)
        for gtxt in gs:
            g = mp.mpf(gtxt)
            Da = lambda h: (1 - g * m.kA(h)) * (1 - g * m.kS(h - 1))
            Ds = lambda h: (1 - g * m.kS(h)) * (1 - g * m.kA(h - 1))
            print(f" q={q} g={gtxt}")
            for kind, k, p in (("A", m.kA, 1), ("S", m.kS, 0)):
                for H in scan_roots(k, g, mp.mpf("0.06"), mp.mpf("3.6")):
                    rad = mp.mpf("0.004")
                    if kind == "A":
                        lit = -4 * _res(lambda h: (2 * h - 1) * psv.k0_psi_s(h) / Da(h), H, rad)
                    else:
                        lit = 4 * _res(lambda h: (2 * h - 1) * psv.k0_psi_a(h) / Ds(h), H, rad)
                    w = 1 / mp.diff(lambda W: RF(m, p, W), H)
                    pred = (-1 if kind == "A" else 1) * 4 * mp.pi * m.norm * m.L * w
                    tag = " (below 1/2)" if H < half else (" (exceptional)" if abs(H - 1) < 1e-8 else "")
                    record(f"{kind} root H={mp.nstr(H, 8)}{tag}: [PSV] weight vs Dirac", lit, pred, tol=mp.mpf(10) ** -7)
                    record(f"   R_F,{p}(H) = g sigma_p L_q and w > 0", RF(m, p, H) * (1 if p == 0 else -1) / (g * m.L), 1, tol=mp.mpf(10) ** -10)
                    RESULTS.append(bool(w > 0))


# ------------------------------------------------------------------ 6
def part_6():
    print("\n6. The mechanism: G_{g,p,11}(w) = -[G^F_{g,p}(w + 1/2) + G^F_{g,1-p}(w - 1/2)] / 2w")
    print(" So the bottom (at w = H - 1/2) and top (at w = H + 1/2) of one Dirac pole have weights in the")
    print(" ratio (H + 1/2) : -(H - 1/2), which is what lets the super-centre integral collapse them.")
    for q in (3, 5):
        m = N2(q)
        for g in (mp.mpf("0.7"), mp.mpf(1)):
            for p in (0, 1):
                GF = lambda pp, W: 1 / (RF(m, pp, W) - g * m.sigma(pp) * m.L)
                for w in (mp.mpf("0.83"), mp.mpc("1.3", "0.4"), mp.mpf("2.71")):
                    G11 = mp.inverse(mp.inverse(m.Dmat(p, w)) - g * m.Mmat(p))[0, 0]
                    record(f"q={q} g={g} p={p} w={mp.nstr(w, 4)}", G11, -(GF(p, w + half) + GF(1 - p, w - half)) / (2 * w), tol=mp.mpf(10) ** -15)


PARTS = {"1": part_1, "2": part_2, "3": part_3, "4": part_4, "5": part_5, "6": part_6}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dps", type=int, default=20)
    parser.add_argument("--parts", default="123456")
    args = parser.parse_args()
    mp.mp.dps = args.dps
    for key in args.parts:
        PARTS[key]()
    failed = RESULTS.count(False)
    print(f"\n{len(RESULTS) - failed}/{len(RESULTS)} checks passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
