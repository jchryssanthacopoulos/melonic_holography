#!/usr/bin/env python3
"""Action of the AdS2 supercharges on normal modes, and the boundary conditions [N2](11), (23).

Charges (paper, Sec. 5), with eps Killing and beta = +i or -i labelling the embedding:

    q(eps, beta):     delta Psi = -M(beta) eps,  delta Phi_b = Psibar eps,
                      delta Phi_psi = -beta Psibar g3 eps,  delta Psibar = 0;
    qbar(ebar, beta): delta Psibar = -ebar Mbar(beta),  delta Phi_b = ebar Psi,
                      delta Phi_psi = conj(beta) ebar g3 Psi,  delta Psi = 0;
    M(beta) = (dslash + mu) Phi_b + beta g3 (dslash - mu) Phi_psi.

All frequencies refer to exp(-i w t). Mode functions are built from the regular hypergeometric
solutions, and boundary data from their connection coefficients.

Usage: python3 scripts/check_bulk_susy_modes.py --dps 20
"""

import argparse
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_robin_identities import N2  # noqa: E402

RESULTS = []
I = mp.mpc(0, 1)
S1 = mp.matrix([[0, 1], [1, 0]])
S3 = mp.matrix([[1, 0], [0, -1]])
G0 = I * S1
G1 = -S3
G3 = G0 * G1
ONE = mp.eye(2)


def record(label, value, reference, tol=None):
    tol = tol if tol is not None else mp.mpf(10) ** (-mp.mp.dps + 9)
    diff = abs(mp.mpc(value) - mp.mpc(reference))
    ok = diff / max(1, abs(mp.mpc(reference))) < tol
    RESULTS.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {mp.nstr(mp.mpc(value), 10)} vs {mp.nstr(mp.mpc(reference), 10)}")


def note(label, value):
    print(f"  [INFO] {label}: {mp.nstr(value, 10)}")


# ---- scalar solutions ------------------------------------------------------------------------


def reg(rho, p, w, th):
    """Regular solution of -f'' + (rho^2 - 1/4) sec^2 f = w^2 f with parity (-1)^p.

    At the right endpoint f -> A zeta^(1/2 - rho) + B zeta^(1/2 + rho).
    """
    half = mp.mpf(1) / 2
    s2 = mp.sin(th) ** 2
    c = mp.cos(th) ** (half + rho)
    if p == 0:
        return c * mp.hyp2f1((half + rho + w) / 2, (half + rho - w) / 2, half, s2)
    return mp.sin(th) * c * mp.hyp2f1((3 * half + rho + w) / 2, (3 * half + rho - w) / 2, 3 * half, s2)


def conn(rho, p, w):
    """(A, B) of reg, from the Gauss connection formula."""
    half = mp.mpf(1) / 2
    if p == 0:
        a, b, c = (half + rho + w) / 2, (half + rho - w) / 2, half
    else:
        a, b, c = (3 * half + rho + w) / 2, (3 * half + rho - w) / 2, 3 * half
    A = mp.gamma(c) * mp.gamma(rho) / (mp.gamma(a) * mp.gamma(b))
    B = mp.gamma(c) * mp.gamma(-rho) / (mp.gamma(c - a) * mp.gamma(c - b))
    return A, B


# ---- Killing data at t=0; e^{-iwt} is restored by the frequency bookkeeping ----------------


def w_A(th):
    return mp.matrix([mp.cos(mp.pi / 4 + th / 2), mp.cos(mp.pi / 4 - th / 2)])


def w_B(th):
    return mp.matrix([-mp.cos(mp.pi / 4 + th / 2), mp.cos(mp.pi / 4 - th / 2)])


def dv(f, th):
    return mp.diff(f, th)


class Model:
    def __init__(self, q):
        self.m = N2(q)
        self.mu, self.d, self.nu = self.m.mu, self.m.delta, self.m.nu

    def dirac_from_v(self, Omega, p):
        """chi(theta) = rot (u, v), with v = Phi_psi parity p at w = Omega - 1/2 ([N2](19))."""
        mu, nu = self.mu, self.nu
        w = Omega - mp.mpf(1) / 2
        v = lambda t: reg(nu, p, w, t)
        u = lambda t: (-dv(v, t) + mu * mp.tan(t) * v(t)) / (w - mu)

        def chi(t):
            a = mp.pi / 4 - t / 2
            uu, vv = u(t), v(t)
            return mp.matrix([mp.cos(a) * uu - mp.sin(a) * vv, mp.sin(a) * uu + mp.cos(a) * vv])

        return chi, u, v

    def T_qbar(self, chi, w_fn, conj_beta):
        """(Phi_b, Phi_psi) from ebar Psi. Radial parts; the factors of zeta cancel."""
        def phib(t):
            return (w_fn(t).T * S1 * chi(t))[0]

        def phip(t):
            return conj_beta * (w_fn(t).T * S1 * G3 * chi(t))[0]

        return phib, phip

    def data(self, phib, phip, w, p):
        """(A_psi, B_psi, A_b, B_b) for a pair of KG solutions with Phi_psi of parity p."""
        nu, d = self.nu, self.d
        t0 = mp.mpf("0.37")
        Cp = phip(t0) / reg(nu, p, w, t0)
        Cb = phib(t0) / reg(d, 1 - p, w, t0)
        Ap, Bp = conn(nu, p, w)
        Ab, Bb = conn(d, 1 - p, w)
        return Cp * Ap, Cp * Bp, Cb * Ab, Cb * Bb

    def bc_residual(self, X, g, p):
        """Residuals of [N2](11) with c = g sigma L."""
        Ap, Bp, Ab, Bb = X
        c = g * self.m.sigma(p) * self.m.L
        return (-2 * self.nu * Bp - c * (2 * self.mu * Ap + 2 * self.d * Bb), Ab - c * Ap)

    def roots(self, branch, g, lo, hi, n=500):
        k = self.m.kA if branch == "A" else self.m.kS
        out = []
        grid = [lo + (hi - lo) * mp.mpf(i) / n for i in range(n + 1)]
        for a, b in zip(grid, grid[1:]):
            fa, fb = k(a) - 1 / g, k(b) - 1 / g
            if fa * fb < 0 and abs(fa) < 20 and abs(fb) < 20:
                out.append(mp.findroot(lambda x: k(x) - 1 / g, (a, b), solver="anderson"))
        return out


def check_lines(md):
    print(" Lines in (A_psi, B_b) of the supercharge images, for a generic Dirac mode")
    mu, d = md.mu, md.d
    half = mp.mpf(1) / 2
    for Omega, p in ((mp.mpf("1.73"), 1), (mp.mpf("2.41"), 0)):
        chi, u, v = md.dirac_from_v(Omega, p)
        # epsilon_A: frequency Omega - 1/2, same parity.
        w = Omega - half
        pb, pp = md.T_qbar(chi, w_A, -I)  # conj_beta = -i, i.e. beta_A = +i
        t0 = mp.mpf("0.41")
        record("A, beta=+i: Phi_b = u", pb(t0), u(t0))
        record("A, beta=+i: Phi_psi = v", pp(t0), v(t0))
        X = md.data(pb, pp, w, p)
        record("A, beta=+i: B_b/A_psi = -(w+mu)/(2 delta)", X[3] / X[0], -(w + mu) / (2 * d))
        # epsilon_B: frequency Omega + 1/2, opposite parity.
        w2 = Omega + half
        for cb, sign, tag in ((I, 1, "beta=-i"), (-I, -1, "beta=+i")):
            pb2, pp2 = md.T_qbar(chi, w_B, cb)
            X2 = md.data(pb2, pp2, w2, 1 - p)
            record(f"B, {tag}: B_b/A_psi = {'+' if sign > 0 else '-'}(w-mu)/(2 delta)", X2[3] / X2[0], sign * (w2 - mu) / (2 * d))


def check_norm(md):
    print(" Norm: chi^T chi = u^2 + v^2 pointwise, so KG norm of the image = 2 w x Dirac norm")
    for Omega, p in ((mp.mpf("1.73"), 1), (mp.mpf("2.41"), 0)):
        chi, u, v = md.dirac_from_v(Omega, p)
        for t0 in (mp.mpf("0.3"), mp.mpf("-1.1")):
            record(f"Omega={Omega} p={p} theta={t0}", (chi(t0).T * chi(t0))[0], u(t0) ** 2 + v(t0) ** 2)


def _psi_from_bosons(md, pb, pp, w, w_fn, beta, t):
    """-M(beta) eps at theta = t, rescaled to psi = Psi/zeta^(1/2): a column vector."""
    mu = md.mu
    sec = 1 / mp.cos(t)
    def piece(f, sgn):
        return G0 * (-I * w) * f(t) + G1 * dv(f, t) + sgn * mu * sec * f(t) * ONE
    return -(piece(pb, 1) + beta * G3 * piece(pp, -1)) * w_fn(t)


def _psibar_from_bosons(md, pb, pp, w, w_fn, conj_beta, t):
    """-ebar Mbar(beta) at theta = t, as a row vector, rescaled the same way."""
    mu = md.mu
    sec = 1 / mp.cos(t)
    def dag(f):
        return -I * S1 * (-I * w) * f(t) - S3 * dv(f, t)
    Mdag = dag(pb) + mu * sec * pb(t) * ONE + conj_beta * (dag(pp) - mu * sec * pp(t) * ONE) * G3
    return -(w_fn(t).T * S1) * (S1 * Mdag * S1)


def check_maps(md):
    print(" Boson -> fermion maps and the raising anticommutator")
    half = mp.mpf(1) / 2
    for w, p in ((mp.mpf("1.63"), 1), (mp.mpf("0.87"), 0), (mp.mpf("2.74"), 1)):
        _check_maps_at(md, w, p)


def _check_maps_at(md, w, p):
    half = mp.mpf(1) / 2
    t0 = mp.mpf("0.41")
    # bottom line at (w, p), from Dirac Omega = w + 1/2, parity p, eps_A with beta = +i
    chi1, _, _ = md.dirac_from_v(w + half, p)
    bot_b, bot_p = md.T_qbar(chi1, w_A, -I)
    # top line at (w, p), from Dirac Omega = w - 1/2, parity 1 - p, eps_B with beta = -i
    chi2, _, _ = md.dirac_from_v(w - half, 1 - p)
    top_b, top_p = md.T_qbar(chi2, w_B, I)
    # q(eps_A, +i) raises to Omega = w + 1/2
    on_top = _psi_from_bosons(md, top_b, top_p, w, w_A, I, t0)
    on_bot = _psi_from_bosons(md, bot_b, bot_p, w, w_A, I, t0)
    record("q(A,+i) kills the top line", abs(on_top[0]) + abs(on_top[1]), 0, tol=mp.mpf(10) ** -10)
    ref = chi1(t0)
    ratio = on_bot[0] / ref[0]
    record("q(A,+i) on the bottom line is proportional to chi", on_bot[1] / ref[1], ratio, tol=mp.mpf(10) ** -10)
    note("   coefficient, against 2 w = " + mp.nstr(2 * w, 6), ratio)
    # qbar(ebar_B, beta_B) also raises; beta_B = -i means conj_beta = +i
    r_top = _psibar_from_bosons(md, top_b, top_p, w, w_B, I, t0)
    record("qbar(B,-i) kills the top line", abs(r_top[0]) + abs(r_top[1]), 0, tol=mp.mpf(10) ** -10)
    r_top_other = _psibar_from_bosons(md, top_b, top_p, w, w_B, -I, t0)
    note("   qbar(B,+i) on the same top line, size (not zero: other embedding)", abs(r_top_other[0]) + abs(r_top_other[1]))
    # anticommutator {q(A,beta_A), qbar(B,beta_B)} on Phi_b: - ebar_B (M - Mbar) eps_A, radial
    def anticomm(pb, pp, beta_A, cb):
        psi = _psi_from_bosons(md, pb, pp, w, w_A, beta_A, t0)
        psib = _psibar_from_bosons(md, pb, pp, w, w_B, cb, t0)
        path1 = (w_B(t0).T * S1 * psi)[0]  # ebar_B Psi
        path2 = (psib * w_A(t0))[0]  # Psibar eps_A
        return path1 - path2
    Eplus = lambda f: -2 * (w_B(t0).T * S1 * (G0 * (-I * w) * f(t0) + G1 * dv(f, t0)) * w_A(t0))[0]
    record("uniform su(1,1|1): {q(A,+i), qbar(B,+i)} Phi_b = E_+ Phi_b, on the bottom line", anticomm(bot_b, bot_p, I, -I), Eplus(bot_b), tol=mp.mpf(10) ** -10)
    record("uniform su(1,1|1): same, on the top line", anticomm(top_b, top_p, I, -I), Eplus(top_b), tol=mp.mpf(10) ** -10)
    record("mixed: {q(A,+i), qbar(B,-i)} kills the top line", anticomm(top_b, top_p, I, I), 0, tol=mp.mpf(10) ** -10)
    v = anticomm(bot_b, bot_p, I, I)
    note("   mixed: {q(A,+i), qbar(B,-i)} on the bottom line", v)
    note("   for comparison, E_+ Phi_b on the bottom line", Eplus(bot_b))
    r_bot_A = _psibar_from_bosons(md, bot_b, bot_p, w, w_A, -I, t0)
    record("qbar(A,+i) kills the bottom line (nilpotency)", abs(r_bot_A[0]) + abs(r_bot_A[1]), 0, tol=mp.mpf(10) ** -10)
    # mode dependence of the mixed anticommutator, in units of the Phi_b part of the top image
    chi3, _, _ = md.dirac_from_v(w + half, p)
    top_of_chi, _ = md.T_qbar(chi3, w_B, I)
    record(f"w={w} p={p}: mixed anticommutator = -4 w (Phi_b of T_B chi)", v / top_of_chi(t0), -4 * w, tol=mp.mpf(10) ** -9)


def check_boundary_form(md):
    """The charges on boundary data, and [N2](11) as a superfield condition.

    Top components on the cylinder: B_minus = 2 delta B_b + mu A_psi,
    B_plus = -2 nu B_psi - mu A_b. Then [N2](11) reads (A_b, B_plus) = c (A_psi, B_minus), and
    with b_F = c a_F it is Phi_+ = c Phi_-, the Amsel-Marolf linear condition on the cylinder.
    Claim: on e^(s tau) with s = -w, q(A,+i) acts as
        a_F = N (s A_psi + B_minus),   b_F = N (s A_b + B_plus).
    q(B,-i), which flips parity, acts with the relative sign needed by c -> -c.
    """
    print(" Boundary form of the charges")
    half = mp.mpf(1) / 2
    mu, d, nu = md.mu, md.d, md.nu
    t0 = mp.mpf("0.41")
    # the two rows of the boundary maps, against the response
    for Omega, pF in ((mp.mpf("1.73"), 1), (mp.mpf("2.41"), 0)):
        RF = 2 * nu * (conn(nu, pF, Omega - half)[1] / conn(nu, pF, Omega - half)[0]) / (Omega - nu)
        chi, _, _ = md.dirac_from_v(Omega, pF)
        pb, pp = md.T_qbar(chi, w_A, -I)
        X = md.data(pb, pp, Omega - half, pF)
        record(f"Omega={Omega}: first row, A_b/A_psi = R_F", X[2] / X[0], RF)
        record(f"Omega={Omega}: first row, B_psi/A_psi = (w-mu) R_F/(2 nu)", X[1] / X[0], (Omega - half - mu) * RF / (2 * nu))
        pb2, pp2 = md.T_qbar(chi, w_B, I)
        X2 = md.data(pb2, pp2, Omega + half, 1 - pF)
        record(f"Omega={Omega}: second row, A_b/A_psi = -R_F", X2[2] / X2[0], -RF)
        record(f"Omega={Omega}: second row, B_psi/A_psi = (w'+mu) R_F/(2 nu)", X2[1] / X2[0], (Omega + half + mu) * RF / (2 * nu))
    # [N2](11) in superfield form, as an algebraic identity
    c, A, B = mp.mpf("0.83"), mp.mpf("1.7"), mp.mpf("-0.6")
    Bpsi_from_11 = -c * (2 * mu * A + 2 * d * B) / (2 * nu)
    Ab = c * A
    Bminus = 2 * d * B + mu * A
    Bplus = -2 * nu * Bpsi_from_11 - mu * Ab
    record("(11) is B_plus = c B_minus", Bplus, c * Bminus)
    for w, p in ((mp.mpf("1.63"), 1), (mp.mpf("2.37"), 0)):
        Ap0, Bp0 = conn(nu, p, w)
        Ab0, Bb0 = conn(d, 1 - p, w)
        # a generic pair: x Phi_psi^reg + y Phi_b^reg
        norm = None
        for x, y in ((mp.mpf(1), mp.mpf("0.37")), (mp.mpf("0.21"), mp.mpf("-1.4"))):
            pp = lambda t, x=x: x * reg(nu, p, w, t)
            pb = lambda t, y=y: y * reg(d, 1 - p, w, t)
            Apsi, Bpsi, Abb, Bbb = x * Ap0, x * Bp0, y * Ab0, y * Bb0
            Bm = 2 * d * Bbb + mu * Apsi
            Bpl = -2 * nu * Bpsi - mu * Abb
            sval = -w
            # q(A,+i): the image Dirac mode has Omega = w + 1/2 and parity p
            out = _psi_from_bosons(md, pb, pp, w, w_A, I, t0)
            chi, u, v = md.dirac_from_v(w + half, p)
            C = out[0] / chi(t0)[0]
            record(f"w={w} p={p}: image is a Dirac mode", out[1], C * chi(t0)[1], tol=mp.mpf(10) ** -9)
            w_D = w  # the image's Phi_psi component has frequency Omega - 1/2 = w
            a0, b0 = conn(nu, p, w_D)
            aF = C * a0
            bF = C * 2 * nu * b0 / (w_D - mu)
            N1 = aF / (sval * Apsi + Bm)
            N2v = bF / (sval * Abb + Bpl)
            if norm is None:
                norm = N1
                print(f"    normalization N = {mp.nstr(N1, 10)}")
            record(f"w={w} p={p}: a_F = N (s A_psi + B_minus)", N1, norm, tol=mp.mpf(10) ** -9)
            record(f"w={w} p={p}: b_F = N (s A_b + B_plus)", N2v, norm, tol=mp.mpf(10) ** -9)


def check_bc(md, g):
    print(f" Boundary conditions at g={g}: images of root Dirac modes")
    half = mp.mpf(1) / 2
    for branch, pF in (("A", 1), ("S", 0)):
        for H in md.roots(branch, g, mp.mpf("1.3"), mp.mpf(5))[:2]:
            chi, u, v = md.dirac_from_v(H, pF)
            pb, pp = md.T_qbar(chi, w_A, -I)
            X = md.data(pb, pp, H - half, pF)
            res = md.bc_residual(X, g, pF)
            record(f"{branch} H={mp.nstr(H, 6)}: eps_A, beta=+i, (11) first", res[0] / X[0], 0, tol=mp.mpf(10) ** -8)
            record(f"{branch} H={mp.nstr(H, 6)}: eps_A, beta=+i, (11) second", res[1] / X[0], 0, tol=mp.mpf(10) ** -8)
            pb2, pp2 = md.T_qbar(chi, w_B, I)
            X2 = md.data(pb2, pp2, H + half, 1 - pF)
            res2 = md.bc_residual(X2, g, 1 - pF)
            record(f"{branch} H={mp.nstr(H, 6)}: eps_B, beta=-i, (11) first", res2[0] / X2[0], 0, tol=mp.mpf(10) ** -8)
            record(f"{branch} H={mp.nstr(H, 6)}: eps_B, beta=-i, (11) second", res2[1] / X2[0], 0, tol=mp.mpf(10) ** -8)
            pb3, pp3 = md.T_qbar(chi, w_B, -I)
            X3 = md.data(pb3, pp3, H + half, 1 - pF)
            res3 = md.bc_residual(X3, g, 1 - pF)
            violated = abs(res3[0] / X3[0]) + abs(res3[1] / X3[0])
            note(f"{branch} H={mp.nstr(H, 6)}: eps_B with the OTHER embedding beta=+i, |violation|", violated)
            RESULTS.append(violated > mp.mpf(10) ** -3)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dps", type=int, default=20)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    for q in (3, 5):
        md = Model(q)
        print(f"\nq={q}")
        check_lines(md)
        check_norm(md)
        check_maps(md)
        check_boundary_form(md)
        for g in (mp.mpf("0.7"), mp.mpf(1)):
            check_bc(md, g)
    failed = RESULTS.count(False)
    print(f"\n{len(RESULTS) - failed}/{len(RESULTS)} checks passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
