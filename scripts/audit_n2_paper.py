#!/usr/bin/env python3
"""Independent audit of research/pdf/n2_dirac_ladders.pdf.

The paper's checks mostly compare closed forms. This script attacks the claims by other routes,
with its own implementations of the kernels and mode functions, direct integrals in place of closed
forms, and tests of statements that were previously derived by hand only. Tests:

  A. Foundations: boundary data against a direct fit; the chain integral in its convergent range.
  B. The reduced kernel by direct integration over y at x = +1 and x = -1, not through A(H), B(H).
  C. Dirac correlators: the spectral sum of free modes against the closed forms, position space.
  D. Closure on the FERMION, which the paper never tested: {q, qbar} against the spinor Lie
     derivative plus an R-rotation; {q, q} = 0; closure on Phi_psi.
  E. Boundary form of qbar (the second line of the paper's (6.10)).
  F. E_+ on bulk modes, and the mismatch in the proof of Prop. 6.1.
  G. Bulycheva's Xi^{AS} from its integral definition (4.24), against (4.39) of that paper.
  H. The generators of that paper with Q != 0, against (A.10)-(A.12). N=1 eigenvalue from the MSW generators.
  I. Poles of Phi^s, Phi^a on ALL chi: the paper says chi > 1 only.
  J. The seed: direct centre integral of the odd part of [PSV](4.9) against the residue formula.
  K. All four entries of [N2](12) from [PSV](4.26), (4.35), which the paper left untested.
  L. The family at g -> 1 on each side of chi = 1 (printed with part I).
  M. Small identities: the flux in superspace, I_a from Thong's weights.
  N. Neutral kernel on (f, -df) in position space; the boson rail; boundary forms of the scalars and
     the Dirac field; the flux and superpotential with fermions; D_psi on the imaginary axis.
  O. Fock space of one multiplet: {Q, Qbar} against L_0 - J_0/2, vacuum energy, the invariant state.
  P. Bulycheva (4.99): printed and corrected factors, summed with the full block, against chi^(2 Delta).

Usage: python3 scripts/audit_n2_paper.py --dps 20 [--parts ABCD...]
"""

import argparse
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

RESULTS = []
I = mp.mpc(0, 1)
half = mp.mpf(1) / 2


def record(label, value, reference, tol=None):
    tol = tol if tol is not None else mp.mpf(10) ** (-mp.mp.dps + 10)
    diff = abs(mp.mpc(value) - mp.mpc(reference))
    ok = diff / max(1, abs(mp.mpc(reference))) < tol
    RESULTS.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}: {mp.nstr(mp.mpc(value), 10)} vs {mp.nstr(mp.mpc(reference), 10)}")
    return ok


def info(label, value):
    print(f"  [INFO] {label}: {mp.nstr(value, 10) if not isinstance(value, str) else value}")


class P:
    """Parameters and fresh definitions, not imported from the other scripts."""

    def __init__(self, q):
        self.q = q
        self.d = mp.mpf(1) / q
        self.nu = 1 - self.d
        self.mu = half - self.d
        self.L = 2 * (q - 1) * mp.power(4, -self.nu)

    def kA(self, H):
        d = self.d
        return (
            -(self.nu * mp.gamma(self.nu) / (mp.pi * d * mp.gamma(d)))
            * mp.gamma(d - H) * mp.gamma(d + H) * (mp.sin(mp.pi * d) - mp.sin(mp.pi * H))
        )

    def kS(self, H):
        d = self.d
        return (
            -(self.nu * mp.gamma(self.nu) / (mp.pi * d * mp.gamma(d)))
            * mp.gamma(d - H) * mp.gamma(d + H) * (mp.sin(mp.pi * d) + mp.sin(mp.pi * H))
        )


def reg(rho, p, w, th):
    s2 = mp.sin(th) ** 2
    c = mp.cos(th) ** (half + rho)
    if p == 0:
        return c * mp.hyp2f1((half + rho + w) / 2, (half + rho - w) / 2, half, s2)
    return mp.sin(th) * c * mp.hyp2f1((3 * half + rho + w) / 2, (3 * half + rho - w) / 2, 3 * half, s2)


def conn(rho, p, w):
    if p == 0:
        a, b, c = (half + rho + w) / 2, (half + rho - w) / 2, half
    else:
        a, b, c = (3 * half + rho + w) / 2, (3 * half + rho - w) / 2, 3 * half
    return (
        mp.gamma(c) * mp.gamma(rho) / (mp.gamma(a) * mp.gamma(b)),
        mp.gamma(c) * mp.gamma(-rho) / (mp.gamma(c - a) * mp.gamma(c - b)),
    )


def resp(rho, p, w):
    A, B = conn(rho, p, w)
    return B / A


def RF(m, p, W):
    return 2 * m.nu * resp(m.nu, p, W - half) / (W - m.nu)


S1 = mp.matrix([[0, 1], [1, 0]])
S3 = mp.matrix([[1, 0], [0, -1]])
ONE = mp.eye(2)
G0 = I * S1
G1 = -S3
G3 = G0 * G1


def wA(t):
    return mp.matrix([mp.cos(mp.pi / 4 + t / 2), mp.cos(mp.pi / 4 - t / 2)])


def wB(t):
    return mp.matrix([-mp.cos(mp.pi / 4 + t / 2), mp.cos(mp.pi / 4 - t / 2)])


def dirac(m, W, p):
    """Radial chi(theta) of the Dirac mode with frequency W and parity p, from (first-order)."""
    w = W - half
    v = lambda t: reg(m.nu, p, w, t)
    u = lambda t: (-mp.diff(v, t) + m.mu * mp.tan(t) * v(t)) / (w - m.mu)

    def chi(t):
        a = mp.pi / 4 - t / 2
        uu, vv = u(t), v(t)
        return mp.matrix([mp.cos(a) * uu - mp.sin(a) * vv, mp.sin(a) * uu + mp.cos(a) * vv])

    return chi


# ------------------------------------------------------------------ A
def part_A():
    print("\nA. Foundations")
    m = P(3)
    for rho, p, w in ((m.nu, 1, mp.mpf("1.3")), (m.d, 0, mp.mpf("0.9"))):
        A, B = conn(rho, p, w)
        # subtract the leading terms; the remainder must be O(zeta^(5/2 - rho))
        for z in (mp.mpf("0.01"), mp.mpf("0.005")):
            th = mp.pi / 2 - z
            val = reg(rho, p, w, th)
            lead = A * z ** (half - rho) + B * z ** (half + rho)
            info(f"rho={mp.nstr(rho, 3)} p={p} zeta={z}: (reg-lead)/zeta^(5/2-rho)", (val - lead) / z ** (5 * half - rho))
        z = mp.mpf("0.003")
        th = mp.pi / 2 - z
        record(f"reg -> A z^(1/2-rho) + B z^(1/2+rho), rho={mp.nstr(rho, 3)}", reg(rho, p, w, th), A * z ** (half - rho) + B * z ** (half + rho), tol=mp.mpf(10) ** -3)
        f = lambda t: reg(rho, p, w, t)
        t0 = mp.mpf("0.4")
        ode = -mp.diff(f, t0, 2) + (rho**2 - mp.mpf(1) / 4) / mp.cos(t0) ** 2 * f(t0) - w**2 * f(t0)
        record("reg solves the ODE", ode, 0, tol=mp.mpf(10) ** -8)
    print(" Chain integral in its convergent range (delta in (1/2, 1)), then continuation is analytic")
    for d in (mp.mpf("0.6"), mp.mpf("0.8")):
        t = mp.mpf("1.7")
        nu = 1 - d
        # written region by region, without the cancellation in z - t
        grid = [0, 1, 3, 10, 30, 100, 300, 800]
        full = [-g_ for g_ in reversed(grid)]
        out = lambda v: (t * mp.exp(v)) ** (1 - d) * (t + t * mp.exp(v)) ** (-d)  # both outer regions, equal by symmetry
        outer = 2 * (mp.quad(out, full) + mp.quad(out, grid))
        inner = lambda v: (t * mp.exp(v)) ** (1 - d) * (t - t * mp.exp(v)) ** (-d)
        mid = -2 * (mp.quad(inner, full) + mp.quad(inner, [0, mp.log(half)]) * 0 + 0) if False else -2 * mp.quad(inner, full[:-1] + [-1, mp.log(half)])
        val = outer + mid
        C = 2 * mp.beta(nu, 2 * d - 1) - mp.beta(nu, nu)
        record(f"delta={d}", val, C * t ** (1 - 2 * d), tol=mp.mpf(10) ** -9)


# ------------------------------------------------------------------ B
def part_B():
    print("\nB. The reduced kernel by direct integration in y")
    for q in (3, 5):
        m = P(q)
        d, nu, mu = m.d, m.nu, m.mu
        # J = 1; rho(t) = rho0 sgn t |t|^(-2 delta), V(y) = v sgn y |y|^(-2 mu); c = v rho0
        b = mp.power(mp.tan(mp.pi * d / 2) / (2 * mp.pi), d)
        C = 2 * mp.beta(nu, 2 * d - 1) - mp.beta(nu, nu)
        rho0 = 2 * mu * b**2 * C
        v = (q - 1) * b ** (q - 2)
        for H in (mp.mpf("0.4") * d, mp.mpf("-0.55") * d):
            for p, ref in ((1, m.kA(H)), (0, m.kS(H))):
                sig = 1 if p == 0 else -1
                for x in (mp.mpf(1), mp.mpf(-1)):
                    def integrand(y):
                        Y = 1 if y > 0 else sig
                        return (
                            rho0 * mp.sign(x + y) * abs(x + y) ** (-2 * d)
                            * v * mp.sign(y) * abs(y) ** (-2 * mu)
                            * abs(y) ** (H - d) * Y
                        )

                    # y = +-exp(v): power laws become exponentials; split at the singular point |y| = 1
                    grid = [0, 1, 3, 10, 30, 100, 300, 800]
                    val = mp.mpf(0)
                    for sgn_y in (1, -1):
                        fn = lambda vv, sgn_y=sgn_y: integrand(sgn_y * mp.exp(vv)) * mp.exp(vv)
                        val += mp.quad(fn, [-g_ for g_ in reversed(grid)]) + mp.quad(fn, grid)
                    Yx = 1 if x > 0 else sig
                    record(f"q={q} H={mp.nstr(H, 3)} p={p} x={x}", val, ref * Yx, tol=mp.mpf(10) ** -8)


# ------------------------------------------------------------------ C
def part_C():
    print("\nC. Dirac correlators: free spectral sum against the closed forms in position space")
    for q in (3, 5):
        m = P(q)
        d = m.d
        cfrak = mp.gamma(d + half) / (mp.sqrt(mp.pi) * mp.gamma(d))

        def spectral(p, t, nmax=200):
            # poles of 1/R_F at H: p=1: delta+1+2n, -(delta+2n); p=0: delta+2n, -(delta+1+2n)
            total = mp.mpf(0)
            for n in range(nmax):
                if p == 1:
                    pos, neg = d + 1 + 2 * n, -(d + 2 * n)
                else:
                    pos, neg = d + 2 * n, -(d + 1 + 2 * n)
                if t < 0:
                    W = pos
                    w = mp.diff(lambda x: 1 / RF(m, p, x), W + mp.mpf(10) ** -12) if False else None
                    res = mp.quad(lambda s: (1 / RF(m, p, W + mp.mpf("0.01") * mp.expjpi(2 * s))) * mp.mpf("0.01") * mp.expjpi(2 * s), [0, 0.5, 1])
                    total += -res * mp.exp(W * t)
                else:
                    W = neg
                    res = mp.quad(lambda s: (1 / RF(m, p, W + mp.mpf("0.01") * mp.expjpi(2 * s))) * mp.mpf("0.01") * mp.expjpi(2 * s), [0, 0.5, 1])
                    total += res * mp.exp(W * t)
                if abs(mp.exp(W * t)) < mp.mpf(10) ** (-mp.mp.dps):
                    break
            return total

        for t in (mp.mpf("0.9"), mp.mpf("-0.9")):
            g1, g0 = spectral(1, t), spectral(0, t)
            same = (g1 + g0) / 2
            cross = (g1 - g0) / 2
            record(f"q={q} t={t}: same-boundary = c sgn t /|2 sinh t/2|^(2 delta)", same, cfrak * mp.sign(t) * abs(2 * mp.sinh(t / 2)) ** (-2 * d), tol=mp.mpf(10) ** -8)
            record(f"q={q} t={t}: cross-boundary = c /(2 cosh t/2)^(2 delta)", cross, cfrak * (2 * mp.cosh(t / 2)) ** (-2 * d), tol=mp.mpf(10) ** -8)


# ------------------------------------------------------------------ D
def boson_to_psi(m, pb, pp, w, wY, beta, t):
    sec = 1 / mp.cos(t)

    def piece(f, sgn):
        return G0 * (-I * w) * f(t) + G1 * mp.diff(f, t) + sgn * m.mu * sec * f(t) * ONE

    return -(piece(pb, 1) + beta * G3 * piece(pp, -1)) * wY(t)


def lie_spinor(chi, W, Xi_t, Xi_th, lam, t):
    """(xi^a nabla_a + (1/4) F_01 gamma^01) on psi, radial, with xi ~ exp(-i lam t)."""
    tan = mp.tan(t)
    sec2 = 1 / mp.cos(t) ** 2
    # Psi = zeta^(1/2) psi, so nabla_theta on Psi, divided by zeta^(1/2): d_theta - (1/2) tan
    dchi = mp.matrix([mp.diff(lambda s: chi(s)[i], t) for i in range(2)])
    ch = chi(t)
    nab_t = (-I * W) * ch - half * tan * (G0 * G1) * ch  # nabla_t = d_t + (1/2) gamma_t^theta d_theta Lam
    nab_th = dchi - half * tan * ch
    # F_t,theta = d_t xi_theta - d_theta xi_t with xi_t = -sec^2 xi^t, xi_theta = sec^2 xi^th
    f_low_t = lambda s: -Xi_t(s) / mp.cos(s) ** 2
    Ftth = (-I * lam) * sec2 * Xi_th(t) - mp.diff(f_low_t, t)
    F01 = Ftth * mp.cos(t) ** 2
    return Xi_t(t) * nab_t + Xi_th(t) * nab_th + mp.mpf(1) / 4 * F01 * (G0 * G1) * ch


def part_D():
    print("\nD. Closure on the fermion, never tested in the paper")
    s = {"A": 1, "B": -1}
    wf = {"A": wA, "B": wB}
    for q in (3, 5):
        m = P(q)
        print(f" q={q}")
        for W, p in ((mp.mpf("1.73"), 1), (mp.mpf("0.61"), 0)):
            chi = dirac(m, W, p)
            for t0 in (mp.mpf("0.37"),):
                for beta, cb, tag in ((I, -I, "beta=+i"), (-I, I, "beta=-i")):
                    # {q(Y), qbar(X)} on Psi: Psi -> (qbar X) -> bosons -> (q Y) -> Psi
                    for Y, X in (("A", "A"), ("B", "B"), ("A", "B"), ("B", "A")):
                        wX, wY = wf[X], wf[Y]
                        wp = W - mp.mpf(s[X]) / 2
                        pb = lambda t: (wX(t).T * S1 * chi(t))[0]
                        pp = lambda t: cb * (wX(t).T * S1 * G3 * chi(t))[0]
                        out = boson_to_psi(m, pb, pp, wp, wY, beta, t0)
                        # xi^a = ebar_X cos(theta) gamma^a eps_Y, radial; exp(-i lam t), lam=(s_Y-s_X)/2
                        lam = mp.mpf(s[Y] - s[X]) / 2
                        Xt = lambda t: (wX(t).T * S1 * G0 * wY(t))[0]
                        Xth = lambda t: (wX(t).T * S1 * G1 * wY(t))[0]
                        lie = lie_spinor(chi, W, Xt, Xth, lam, t0)
                        # expected: out = -2 Lie + alpha chi, alpha = R-rotation
                        ch = chi(t0)
                        rem = out + 2 * lie
                        alpha0 = rem[0] / ch[0]
                        alpha1 = rem[1] / ch[1]
                        expected = {"AA": 1, "BB": -1}.get(Y + X, 0)
                        if tag == "beta=+i":
                            info(f"   (Y,X)=({Y},{X}) sizes: |out|, |L_xi Psi|, |chi|", f"{mp.nstr(abs(out[0]) + abs(out[1]), 5)}, {mp.nstr(abs(lie[0]) + abs(lie[1]), 5)}, {mp.nstr(abs(ch[0]) + abs(ch[1]), 5)}")
                        ok0 = record(f"W={W} p={p} {tag} (Y,X)=({Y},{X}): {{q,qbar}} Psi = -2 L_xi Psi + alpha Psi, comp 0, alpha", alpha0, expected, tol=mp.mpf(10) ** -7)
                        record(f"W={W} p={p} {tag} (Y,X)=({Y},{X}): comp 1, alpha", alpha1, expected, tol=mp.mpf(10) ** -7)
                        _ = ok0
                    # {q(Y), q(X)} on Psibar -> Psi, must vanish in a uniform embedding
                    for Y, X in (("A", "A"), ("B", "B"), ("A", "B")):
                        total = mp.matrix([0, 0])
                        for a, c in ((Y, X), (X, Y)):
                            wc, wa = wf[c], wf[a]
                            wp = -W + mp.mpf(s[c]) / 2
                            rowchi = lambda t: chi(t).T * S1
                            pb = lambda t, wc=wc: (rowchi(t) * wc(t))[0]
                            pp = lambda t, wc=wc: -beta * (rowchi(t) * G3 * wc(t))[0]
                            total += boson_to_psi(m, pb, pp, wp, wa, beta, t0)
                        record(f"W={W} p={p} {tag}: {{q({Y}),q({X})}} = 0", abs(total[0]) + abs(total[1]), 0, tol=mp.mpf(10) ** -7)
            # the mixed embedding, for information
            total = mp.matrix([0, 0])
            t0 = mp.mpf("0.37")
            for a, c, bc, ba in (("A", "B", -I, I), ("B", "A", I, -I)):
                wc, wa = wf[c], wf[a]
                wp = -W + mp.mpf(s[c]) / 2
                rowchi = lambda t: chi(t).T * S1
                pb = lambda t, wc=wc: (rowchi(t) * wc(t))[0]
                pp = lambda t, wc=wc, bc=bc: -bc * (rowchi(t) * G3 * wc(t))[0]
                total += boson_to_psi(m, pb, pp, wp, wa, ba, t0)
            info(f"W={W} p={p}: mixed {{q(A,+i), q(B,-i)}} size (nonzero means a new bosonic map)", abs(total[0]) + abs(total[1]))
        print(" Closure on Phi_psi, on a generic pair of scalars, uniform embedding")
        w = mp.mpf("1.63")
        p = 1
        t0 = mp.mpf("0.37")
        for x, y in ((mp.mpf(1), mp.mpf("0.37")), (mp.mpf("0.2"), mp.mpf("-1.1"))):
            pp = lambda t, x=x: x * reg(m.nu, p, w, t)
            pb = lambda t, y=y: y * reg(m.d, 1 - p, w, t)
            for beta, cb in ((I, -I), (-I, I)):
                psi = boson_to_psi(m, pb, pp, w, wA, beta, t0)
                sec = 1 / mp.cos(t0)

                def dag(f):
                    return -I * S1 * (-I * w) * f(t0) - S3 * mp.diff(f, t0)

                Mdag = dag(pb) + m.mu * sec * pb(t0) * ONE + cb * (dag(pp) - m.mu * sec * pp(t0) * ONE) * G3
                row = -(wB(t0).T * S1) * (S1 * Mdag * S1)
                path1 = cb * (wB(t0).T * S1 * G3 * psi)[0]
                path2 = -beta * (row * G3 * wA(t0))[0]
                Xt, Xth = (wB(t0).T * S1 * G0 * wA(t0))[0], (wB(t0).T * S1 * G1 * wA(t0))[0]
                Eplus = -2 * (Xt * (-I * w) * pp(t0) + Xth * mp.diff(pp, t0))
                record(f"x={mp.nstr(x, 3)} y={mp.nstr(y, 3)} beta={mp.nstr(beta, 2)}: {{q_A,qbar_B}} Phi_psi = E_+ Phi_psi", path1 - path2, Eplus, tol=mp.mpf(10) ** -8)


# ------------------------------------------------------------------ E
def part_E():
    print("\nE. Boundary form of qbar_A: delta abar = s A - B, on modes s = -w")
    for q in (3, 5):
        m = P(q)
        mu, d, nu = m.mu, m.d, m.nu
        t0 = mp.mpf("0.41")
        for w, p in ((mp.mpf("1.63"), 1), (mp.mpf("2.37"), 0)):
            Ap0, Bp0 = conn(nu, p, w)
            Ab0, Bb0 = conn(d, 1 - p, w)
            norm = None
            for x, y in ((mp.mpf(1), mp.mpf("0.37")), (mp.mpf("0.21"), mp.mpf("-1.4"))):
                pp = lambda t, x=x: x * reg(nu, p, w, t)
                pb = lambda t, y=y: y * reg(d, 1 - p, w, t)
                sec = 1 / mp.cos(t0)

                def dag(f):
                    return -I * S1 * (-I * w) * f(t0) - S3 * mp.diff(f, t0)

                # beta = +i, so conj(beta) = -i in Mbar
                Mdag = dag(pb) + mu * sec * pb(t0) * ONE + (-I) * (dag(pp) - mu * sec * pp(t0) * ONE) * G3
                row = -(wA(t0).T * S1) * (S1 * Mdag * S1)
                # Psibar = Psi^dagger sigma_1, so psi^dagger = row sigma_1; it is the conjugate of a mode
                # of frequency W' = -(w - 1/2), with real radial function
                dag_psi = row * S1
                Wp = -(w - half)
                chi = dirac(m, Wp, p)
                ch = chi(t0)
                C = dag_psi[0] / ch[0]
                record(f"w={w} p={p}: image is conjugate of a Dirac mode", dag_psi[1], C * ch[1], tol=mp.mpf(10) ** -8)
                wD = Wp - half
                a0, b0 = conn(nu, p, wD)
                abar = C * a0
                bbar = C * 2 * nu * b0 / (wD - mu)
                Apsi, Bpsi, Abb, Bbb = x * Ap0, x * Bp0, y * Ab0, y * Bb0
                Bm = 2 * d * Bbb + mu * Apsi
                Bpl = -2 * nu * Bpsi - mu * Abb
                s = -w
                N1 = abar / (s * Apsi - Bm)
                N2 = bbar / (s * Abb - Bpl)
                if norm is None:
                    norm = N1
                    info(f"w={w} p={p}: normalization", norm)
                record(f"w={w} p={p}: abar = N (s A_psi - B_minus)", N1, norm, tol=mp.mpf(10) ** -8)
                record(f"w={w} p={p}: bbar = N (s A_b - B_plus)", N2, norm, tol=mp.mpf(10) ** -8)


# ------------------------------------------------------------------ F
def part_F():
    print("\nF. E_+ on bulk modes, and the mismatch in Prop. 6.1")
    m = P(3)
    for rho, p, w in ((m.nu, 1, mp.mpf("1.63")), (m.d, 0, mp.mpf("1.63")), (m.nu, 0, mp.mpf("0.9"))):
        f = lambda t: reg(rho, p, w, t)
        Ef = lambda t: w * mp.sin(t) * f(t) - mp.cos(t) * mp.diff(f, t)
        t0 = mp.mpf("0.4")
        C = Ef(t0) / reg(rho, 1 - p, w + 1, t0)
        record(f"rho={mp.nstr(rho, 3)} p={p}: E_+ f solves KG at w+1, parity flipped", Ef(mp.mpf("0.9")), C * reg(rho, 1 - p, w + 1, mp.mpf("0.9")), tol=mp.mpf(10) ** -8)
        A0, B0 = conn(rho, p, w)
        A1, B1 = conn(rho, 1 - p, w + 1)
        kA_, kB_ = C * A1 / A0, C * B1 / B0
        record("A' = kappa (w + Delta_-) A", kA_ / (w + half - rho), kB_ / (w + half + rho), tol=mp.mpf(10) ** -8)
        info("   kappa", kA_ / (w + half - rho))
    print(" Mismatch: a BC mode, transformed, against -2 g w (M_p)_21 (X0)_1")
    g = mp.mpf("0.7")
    for p, k in ((1, m.kA), (0, m.kS)):
        H = mp.findroot(lambda x: k(x) - 1 / g, mp.mpf("3.1") if p == 1 else mp.mpf("2.67"))
        w = H - half
        sigma = 1 if p == 0 else -1
        c = g * sigma * m.L
        Apsi = mp.mpf(1)
        Bb = -(w + m.mu) / (2 * m.d) * Apsi
        Ab = c * Apsi
        Bpsi = (w - m.mu) * Ab / (2 * m.nu)
        # transform with weight Delta: X' = -(w + Delta) X
        Ap_, Bb_ = -(w + -m.mu) * Apsi, -(w + 1 - m.mu) * Bb
        Ab_, Bp_ = -(w + m.mu) * Ab, -(w + 1 + m.mu) * Bpsi
        cp = -c
        res2 = 2 * m.d * Ab_ - cp * 2 * m.d * Ap_
        record(f"p={p}: second component of j0' - g M_(1-p) X'", res2, -2 * w * c * 2 * m.d * Apsi, tol=mp.mpf(10) ** -8)


# ------------------------------------------------------------------ G
def part_G():
    print("\nG. Bulycheva's Xi^AS from the integral (4.24), against (4.39)")

    def Fb(h, z):
        return mp.gamma(h) ** 2 / mp.gamma(2 * h) * z**h * mp.hyp2f1(h, h, 2 * h, z)

    def phi(h, z):
        return Fb(h, z) - Fb(h + 1, z)

    for h in (mp.mpf("0.3"), mp.mpf("0.7")):
        for chi_ in (mp.mpf("0.25"), mp.mpf("0.6")):
            def Q(y, dist_chi):
                # h chi^h (-sgn(1-y)) |1-y|^h |chi-y|^(-h), with |chi-y| passed exactly
                return h * chi_**h * (-mp.sign(1 - y)) * abs(1 - y) ** h * dist_chi ** (-h)

            def f(y):
                return (
                    Q(y, abs(chi_ - y)) * abs(y) ** (-h) * (1 / y + 1 / (chi_ - y) - 1 / (1 - y))
                )

            # principal values, folded by hand so that the small distances are exact, not differences
            def fold0(e):
                plus = Q(e, chi_ - e) * e ** (-h) * (1 / (chi_ - e) - 1 / (1 - e))
                minus_ = Q(-e, chi_ + e) * e ** (-h) * (1 / (chi_ + e) - 1 / (1 + e))
                diff = (Q(e, chi_ - e) - Q(-e, chi_ + e)) * e ** (-h) / e
                return plus + minus_ + diff

            def foldc(e):
                yp, ym = chi_ + e, chi_ - e
                gp = h * chi_**h * (-mp.sign(1 - yp)) * abs(1 - yp) ** h * yp ** (-h)
                gm = h * chi_**h * (-mp.sign(1 - ym)) * abs(1 - ym) ** h * ym ** (-h)
                reg_ = e ** (-h) * (gp * (1 / yp - 1 / (1 - yp)) + gm * (1 / ym - 1 / (1 - ym)))
                sing = e ** (-h) * (gm - gp) / e
                return reg_ + sing

            a = chi_ / 2
            b = min(chi_ / 2, (1 - chi_) / 2)
            with mp.workdps(45):
                pv0 = mp.quad(fold0, [0, a / 4, a], maxdegree=8)
                pvc = mp.quad(foldc, [0, b / 4, b], maxdegree=8)
                rest = mp.quad(f, [-mp.inf, -4, -1, -a], maxdegree=8)
                if chi_ - b > a:
                    rest += mp.quad(f, [a, chi_ - b], maxdegree=8)
                rest += mp.quad(f, [chi_ + b, (chi_ + b + 1) / 2, 1], maxdegree=8)
                rest += mp.quad(f, [1, 2, 8, mp.inf], maxdegree=8)
            val = pv0 + pvc + rest
            sec = 1 / mp.cos(mp.pi * h)
            minus = h * ((1 + sec) * phi(h, chi_) - (1 - sec) * phi(-h, chi_))
            plus = h * ((1 + sec) * phi(h, chi_) + (1 - sec) * phi(-h, chi_))
            record(f"h={h} chi={chi_}: (4.24) = corrected (4.39), minus sign", val, minus, tol=mp.mpf(10) ** -6)
            info("   printed (4.39), plus sign, would give", plus)


# ------------------------------------------------------------------ H
def part_H():
    print("\nH. Generators with Q != 0 against (A.10)-(A.12), and N=1")
    import sympy as sp

    from check_supercasimir import G as GG
    from check_supercasimir import angle, casimir, compose, generators, lift, odd, size, d_even, d_odd, total, POINT

    _ = (GG, angle)
    D = sp.Rational(1, 6)
    Q = sp.Rational(2, 7)
    tau = sp.Symbol("tau1", real=True)
    for fac, label in ((sp.Integer(1), "Q/2 -> Q"), (sp.Rational(1, 2), "printed Q/2")):
        g = generators(1, D, Q, fac)
        cas = casimir(g, "right")
        f = lift(tau ** sp.Rational(3, 5))
        out = cas(f)
        tar = (D**2 - Q**2 / 4) * f
        print(f"    {label}: C f(tau) vs (D^2 - Q^2/4) f: residual {sp.N(size(out - tar), 6)}")
        th = odd("t1") * lift(tau ** sp.Rational(3, 5))
        tb = odd("tb1") * lift(tau ** sp.Rational(3, 5))
        print(f"    {label}: C (theta f) coefficient ratio", sp.N((cas(th).terms.get(("t1",), 0) / (tau ** sp.Rational(3, 5))).subs(POINT), 8), " (A.11) gives", sp.N(D**2 - Q**2 / 4 + Q / 4, 8))
        print(f"    {label}: C (thetabar f) coefficient ratio", sp.N((cas(tb).terms.get(("tb1",), 0) / (tau ** sp.Rational(3, 5))).subs(POINT), 8), " (A.12) gives", sp.N(D**2 - Q**2 / 4 - Q / 4, 8))

    # alternative repair: keep Q/2 in (A.3), (A.5), (A.7) and put Q/2 in J0
    print(" N=1: the MSW generators on f(t12 - theta1 theta2)")
    half_ = sp.Rational(1, 2)
    Dl = sp.Rational(1, 6)

    def gens1(n):
        t = sp.Symbol(f"tau{n}", real=True)
        th = odd(f"t{n}")
        dt = lambda f: d_even(f, t)
        dth = lambda f: d_odd(f, f"t{n}")
        return {
            "L0": lambda f: -t * dt(f) - half_ * (th * dth(f)) - Dl * f,
            "Lm": lambda f: -dt(f),
            "Lp": lambda f: -t**2 * dt(f) - t * (th * dth(f)) - 2 * Dl * t * f,
            "Gm": lambda f: dth(f) - th * dt(f),
            "Gp": lambda f: t * dth(f) - t * (th * dt(f)) - 2 * Dl * (th * f),
        }

    gs = (gens1(1), gens1(2))
    tot = {k: total(k, gs) for k in ("L0", "Lm", "Lp", "Gm", "Gp")}
    t1, t2 = sp.Symbol("tau1", real=True), sp.Symbol("tau2", real=True)
    x = lift(t1 - t2) - (odd("t1") * odd("t2"))
    for h in (sp.Rational(7, 5), sp.Rational(9, 4)):
        f = compose(lambda u: u ** (h - 2 * Dl), x)
        Cf = (
            tot["L0"](tot["L0"](f))
            - half_ * (tot["Lm"](tot["Lp"](f)) + tot["Lp"](tot["Lm"](f)))
            - sp.Rational(1, 4) * (tot["Gm"](tot["Gp"](f)) - tot["Gp"](tot["Gm"](f)))
        )
        res = size(Cf - h * (h - half_) * f)
        ok = res < 1e-20
        RESULTS.append(bool(ok))
        print(f"  [{'PASS' if ok else 'FAIL'}] h={h}: C f = h (h - 1/2) f, residual {sp.N(res, 4)}")
        print(f"    total G_(-1/2) on f: size {sp.N(size(tot['Gm'](f)), 4)}, total L_(-1): {sp.N(size(tot['Lm'](f)), 4)}")


# ------------------------------------------------------------------ I, L
def part_I():
    print("\nI. Are the poles of Phi^s, Phi^a confined to chi > 1?")

    def phi_s(h, chi):
        if chi > 1:
            x = (chi - 2) / chi
            return mp.gamma(half - h / 2) * mp.gamma(h / 2) / mp.sqrt(mp.pi) * mp.hyp2f1(h / 2, half - h / 2, half, x**2)
        A = half * mp.tan(mp.pi * h) * mp.cot(mp.pi * h / 2)
        B = -half * mp.tan(mp.pi * h) * mp.tan(mp.pi * h / 2)
        F1 = chi**h * mp.gamma(h) ** 2 / mp.gamma(2 * h) * mp.hyp2f1(h, h, 2 * h, chi)
        F2 = chi ** (1 - h) * mp.gamma(1 - h) ** 2 / mp.gamma(2 - 2 * h) * mp.hyp2f1(1 - h, 1 - h, 2 - 2 * h, chi)
        return A * F1 + B * F2

    rad = mp.mpf("0.02")
    for chi_ in (mp.mpf("0.5"), mp.mpf("0.2"), mp.mpf("1.6")):
        res = mp.quad(lambda t: phi_s(1 + rad * mp.expjpi(2 * t), chi_) * rad * mp.expjpi(2 * t), [0, 0.25, 0.5, 0.75, 1])
        info(f"Res_(h=1) Phi^s at chi={chi_}", res)
    print("    So the eigenfunction has the pole on ALL chi > 0; only the F_1 representation on 0<chi<1 has none.")

    print("\nL. The family at g -> 1, on each side of chi = 1 (p=1, q=3)")
    from check_literature_inputs import PSV  # transcription of (4.10), checked elsewhere

    m = P(3)
    psv = PSV(3)
    for g in (mp.mpf("0.99"), mp.mpf("0.999"), mp.mpf("0.9999")):
        H = mp.findroot(lambda x: m.kA(x) - 1 / g, mp.mpf(1) + (1 - g) * mp.mpf("0.1"), tol=mp.mpf(10) ** -15)
        D = lambda h: (1 - g * m.kA(h)) * (1 - g * m.kS(h - 1))
        two_mu = lambda h: (2 * h - 1) / (mp.pi * mp.tan(mp.pi * h / 2))
        rad = abs(H - 1) / 4
        # coefficient of the block on each side: residue of 2 mu k0 / D times the function there
        Rk = mp.quad(lambda t: (psv.k0_psi_s(H + rad * mp.expjpi(2 * t)) / D(H + rad * mp.expjpi(2 * t))) * rad * mp.expjpi(2 * t), [0, 0.25, 0.5, 0.75, 1])
        inside = -two_mu(H) * Rk * mp.gamma(H) ** 2 / mp.gamma(2 * H) * mp.mpf("0.5") ** H * mp.hyp2f1(H, H, 2 * H, mp.mpf("0.5"))
        outside = -two_mu(H) * Rk * phi_s(H, mp.mpf("1.6"))
        info(f"g={g}, H={mp.nstr(H, 8)}: F contribution at chi=0.5 / at chi=1.6 / -Res[k0/D]", f"{mp.nstr(inside, 6)} / {mp.nstr(outside, 8)} / {mp.nstr(-Rk, 8)}")
    print("    The last two columns should agree as g -> 1, and the first should vanish.")


# ------------------------------------------------------------------ J
def part_J():
    print("\nJ. The odd seed: a direct centre integral of [PSV](4.9) against the residue formula")
    from check_literature_inputs import PSV

    for q in (3, 5):
        m = P(q)
        psv = PSV(q)
        d = m.d

        def Fs(chi):
            # part of (4.9) symmetric under chi -> chi/(chi-1)
            return half * (mp.sign(chi / (1 - chi)) * abs(chi / (1 - chi)) ** d - mp.sign(chi) * abs(chi) ** d)

        for tau in (mp.mpf("0.6"), mp.mpf("1.4")):
            s, c = mp.sinh(tau / 2), mp.cosh(tau / 2)
            f = lambda u: Fs(1 / (u**2 - s**2))
            direct = 2 * (mp.quad(f, [0, s]) + mp.quad(f, [s, c]) + mp.quad(f, [c, 2 * c, 10 * c, 100 * c, mp.inf]))
            total = mp.mpf(0)
            for n in range(0, 40):
                h0 = 1 + d + 2 * n
                rad = mp.mpf("0.05")
                res = mp.quad(lambda t: psv.k0_psi_s(h0 + rad * mp.expjpi(2 * t)) * rad * mp.expjpi(2 * t), [0, 0.25, 0.5, 0.75, 1])
                total += -2 * res * mp.exp(-(h0 - half) * tau)
            record(f"q={q} tau={tau}: P[F_0, odd] = -2 sum Res k0^(psi,s) e^(-w tau)", direct, total, tol=mp.mpf(10) ** -5)


# ------------------------------------------------------------------ K
def part_K():
    print("\nK. All four entries of [N2](12), from [PSV](4.26), (4.35)")
    from check_literature_inputs import PSV

    def res(f, h0, rad=mp.mpf("0.004")):
        return mp.quad(lambda t: f(h0 + rad * mp.expjpi(2 * t)) * rad * mp.expjpi(2 * t), [0, 0.25, 0.5, 0.75, 1])

    for q in (3, 5):
        m = P(q)
        psv = PSV(q)
        d, nu, mu = m.d, m.nu, m.mu
        b = mp.power(mp.tan(mp.pi * d / 2) / (2 * mp.pi), d)
        bb = d * b
        ell = 2 * mp.power(4, -nu) / b ** (q - 2)
        Sdiag = (mp.mpf(1), 1 / (2 * d))
        bs = (b, bb)

        def k0b_s(h):  # (4.12)
            return (mp.pi * mp.cos(mp.pi / (2 * q)) ** 2 * mp.gamma(-d) ** 2) / ((mp.sin(mp.pi * h) + mp.sin(mp.pi * d)) * mp.gamma(-h - d + 1) * mp.gamma(h - d))

        def k0b_a(h):  # (4.13)
            return (mp.pi * mp.cos(mp.pi / (2 * q)) ** 2 * mp.gamma(-d) ** 2) / ((mp.sin(mp.pi * h) - mp.sin(mp.pi * d)) * mp.gamma(-h - d + 1) * mp.gamma(h - d))

        for g in (mp.mpf("0.9"), mp.mpf(1)):
            for p in (1, 0):
                sigma = 1 if p == 0 else -1
                mat = psv.ka_matrix if p == 1 else psv.ks_matrix
                k0p = psv.k0_psi_s if p == 1 else psv.k0_psi_a
                k0b = k0b_s if p == 1 else k0b_a
                Dg = (lambda h: (1 - g * m.kA(h)) * (1 - g * m.kS(h - 1))) if p == 1 else (lambda h: (1 - g * m.kS(h)) * (1 - g * m.kA(h - 1)))
                sign = -2 if p == 1 else 2

                def Gres(wr):
                    def Gmat(w):
                        Dd = mp.matrix([[-1 / (2 * nu * (conn(nu, p, w)[1] / conn(nu, p, w)[0])), 0], [0, (conn(d, 1 - p, w)[1] / conn(d, 1 - p, w)[0]) / (2 * d)]])
                        M = sigma * m.L * mp.matrix([[2 * mu, 2 * d], [2 * d, 0]])
                        return mp.inverse(mp.inverse(Dd) - g * M)

                    return [[-res(lambda w: Gmat(w)[i, j], wr) for j in range(2)] for i in range(2)]

                # every root with omega < 4: bottom components of this sector's branch, and top components
                # of the other branch at h = H + 1, which are the two kinds of zero of Dg. Includes the
                # exceptional roots h = 1, 2 at g = 1.
                branch1 = m.kA if p == 1 else m.kS
                branch2 = m.kS if p == 1 else m.kA

                def scan_roots(k, lo, hi):
                    out = []
                    x = lo + mp.sqrt(2) * mp.mpf(10) ** -4
                    while x < hi:
                        fa, fb = k(x) - 1 / g, k(x + mp.mpf("0.01")) - 1 / g
                        if fa * fb < 0 and abs(fa) < 30 and abs(fb) < 30:
                            out.append(mp.findroot(lambda y: k(y) - 1 / g, (x, x + mp.mpf("0.01")), solver="anderson"))
                        x += mp.mpf("0.01")
                    return out

                roots_h = scan_roots(branch1, mp.mpf("0.6"), mp.mpf("4.5")) + [H + 1 for H in scan_roots(branch2, mp.mpf("0.6"), mp.mpf("3.5"))]
                for h0 in roots_h:
                    wr = h0 - half
                    W = Gres(wr)
                    entries = {
                        (0, 0): lambda h: k0p(h) / Dg(h),
                        (1, 0): lambda h: g * mat(h)[1, 0] * k0p(h) / Dg(h),
                        (0, 1): lambda h: g * mat(h)[0, 1] * k0b(h) / Dg(h),
                        (1, 1): lambda h: (1 - g * mat(h)[0, 0]) * k0b(h) / Dg(h),
                    }
                    vals = {}
                    for (i, j), fn in entries.items():
                        vals[(i, j)] = sign * res(fn, h0)
                        # the matrix (4.25) acts on coefficients of un-normalized bilocals. So the normalized
                        # correlator of column j carries b_j^2 in both rows; this is forced, see below.
                        lhs = 2 * Sdiag[i] * Sdiag[j] * bs[j] * bs[j] * vals[(i, j)]
                        rhs = sigma * ell * W[i][j]
                        record(f"q={q} g={g} p={p} h={mp.nstr(h0, 6)} entry {i + 1}{j + 1}", lhs, rhs, tol=mp.mpf(10) ** -6)
                    # Literal reading: prefactor G_x(12) G_y(34), amplitude b_i b_j. Symmetry F_21 = F_12 decides.
                    lit21 = bs[1] * bs[0] * vals[(1, 0)]
                    lit12 = bs[0] * bs[1] * vals[(0, 1)]
                    col21 = bs[0] * bs[0] * vals[(1, 0)]
                    col12 = bs[1] * bs[1] * vals[(0, 1)]
                    info(f"   F_21/F_12 at this root: literal prefactor {mp.nstr(mp.re(lit21 / lit12), 6)} (delta^2 = {mp.nstr(d**2, 6)}); column prefactor", mp.re(col21 / col12))
                    RESULTS.append(abs(col21 / col12 - 1) < mp.mpf(10) ** -6)


def part_M():
    print("\nM. Small identities: the flux in superspace; I_a from Thong's weights")
    m = P(3)
    mu, d, nu = m.mu, m.d, m.nu
    # delta S_free = int d^2 theta Phi_+ delta Phi_- reproduces j0 . delta X, for arbitrary data
    A, B, Ab, Bp, dA, dB = [mp.mpf(x) for x in ("0.7", "-1.3", "0.4", "2.1", "0.9", "-0.2")]
    Bm = 2 * d * B + mu * A
    Bpl = -2 * nu * Bp - mu * Ab
    dBm = 2 * d * dB + mu * dA
    record("A_b dB_- + B_+ dA = j0 . dX", Ab * dBm + Bpl * dA, (-2 * nu * Bp) * dA + (2 * d * Ab) * dB)
    c = mp.mpf("0.83")
    # Euler-Lagrange of S = free - (c/2) int d^2 theta Phi_-^2 gives both components of (11)
    coeff_dB = 2 * d * (Ab - c * A)
    coeff_dA = mu * (Ab - c * A) + Bpl - c * Bm
    record("coefficient of dB_b is 2 delta (A_b - c A)", coeff_dB, 2 * d * Ab - c * 2 * d * A)
    record("coefficient of dA_psi is -2 nu B_psi - c (2 mu A + 2 delta B_b)", coeff_dA, -2 * nu * Bp - c * (2 * mu * A + 2 * d * B))
    # Thong's weights with tan(pi h/2) -> tan(pi (h-1)/2) give the even-sector I_a
    for h in (mp.mpf("1.3"), mp.mpf("2.4")):
        w = h - half
        Ta = mp.tan(mp.pi * (h - 1) / 2)
        val = mp.cos(mp.pi * w) + Ta * mp.sin(mp.pi * w)
        record(f"h={h}: cos(pi w) + tan(pi(h-1)/2) sin(pi w) = cot(pi h/2)", val, mp.cot(mp.pi * h / 2))
        T = mp.tan(mp.pi * h / 2)
        record(f"h={h}: cos(pi w) + tan(pi h/2) sin(pi w) = tan(pi h/2), [T](C.8)", mp.cos(mp.pi * w) + T * mp.sin(mp.pi * w), T)


# ------------------------------------------------------------------ E2
def part_E2():
    print("\nE2. Fermion -> boson: delta A and delta B on a regular Dirac mode at generic frequency")
    print(" A mode with free frequency carries generic boundary data, so this is the generic statement.")
    for q in (3, 5):
        m = P(q)
        mu, d, nu = m.mu, m.d, m.nu
        t0, t1 = mp.mpf("0.41"), mp.mpf("0.93")
        for Om, p in ((mp.mpf("2.13"), 1), (mp.mpf("2.87"), 0), (mp.mpf("0.79"), 1)):
            w = Om - half
            chi = dirac(m, Om, p)
            # beta = +i: qbar_A gives (Phi_b, Phi_psi) = (ebar psi, conj(beta) ebar g3 psi), radial
            Pb = lambda t: (wA(t).T * S1 * chi(t))[0]
            Pp = lambda t: -I * (wA(t).T * S1 * G3 * chi(t))[0]
            # q_A on the conjugate mode: (psibar eps, -beta psibar g3 eps), radial
            Pb2 = lambda t: (chi(t).T * S1 * wA(t))[0]
            Pp2 = lambda t: -I * (chi(t).T * S1 * G3 * wA(t))[0]
            record(f"q={q} W={Om} p={p}: q_A and qbar_A give the same radial scalars", Pp2(t0), Pp(t0), tol=mp.mpf(10) ** -10)
            record(f"q={q} W={Om} p={p}: same, Phi_b", Pb2(t0), Pb(t0), tol=mp.mpf(10) ** -10)
            Cb = Pb(t0) / reg(d, 1 - p, w, t0)
            Cp = Pp(t0) / reg(nu, p, w, t0)
            record(f"q={q} W={Om} p={p}: image solves KG, Phi_b", Pb(t1), Cb * reg(d, 1 - p, w, t1), tol=mp.mpf(10) ** -8)
            record(f"q={q} W={Om} p={p}: image solves KG, Phi_psi", Pp(t1), Cp * reg(nu, p, w, t1), tol=mp.mpf(10) ** -8)
            Ap0, Bp0 = conn(nu, p, w)
            Ab0, Bb0 = conn(d, 1 - p, w)
            aF, bF = Ap0, 2 * nu * Bp0 / (w - mu)
            Apsi, Bpsi, Ab, Bb = Cp * Ap0, Cp * Bp0, Cb * Ab0, Cb * Bb0
            Bm, Bpl = 2 * d * Bb + mu * Apsi, -2 * nu * Bpsi - mu * Ab
            record(f"q={q} W={Om} p={p}: delta A_psi = a_F", Apsi, aF, tol=mp.mpf(10) ** -8)
            record(f"q={q} W={Om} p={p}: delta A_b = b_F", Ab, bF, tol=mp.mpf(10) ** -8)
            record(f"q={q} W={Om} p={p}: delta B_- = -w a_F, i.e. d_tau a (qbar) or -d_tau abar (q)", Bm, -w * aF, tol=mp.mpf(10) ** -8)
            record(f"q={q} W={Om} p={p}: delta B_+ = -w b_F", Bpl, -w * bF, tol=mp.mpf(10) ** -8)
    print(" The graded form of (6.10) on (A, B, a, ahat), with d_tau -> s, closes only for ahat = -abar")
    s_ = mp.mpf("0.73")

    def reps(sign):
        # rows: (q phi_i) = sum M_ij phi_j, fields ordered A, B, a, ahat
        Mq = mp.matrix([[0, 0, 0, 1], [0, 0, 0, -s_], [s_, 1, 0, 0], [0, 0, 0, 0]])
        Mb = mp.matrix([[0, 0, 1, 0], [0, 0, s_, 0], [0, 0, 0, 0], [sign * s_, -sign, 0, 0]])
        return Mq, Mb

    Mq, Mb = reps(1)
    An = Mq * Mb + Mb * Mq
    record("q^2 = 0", mp.norm(Mq * Mq), 0)
    record("qbar^2 = 0", mp.norm(Mb * Mb), 0)
    record("{q, qbar} = 2 d_tau on all four fields", mp.norm(An - 2 * s_ * mp.eye(4)), 0)
    Mq, Mb = reps(-1)
    info("with the wrong sign, |{q,qbar} - 2 s| (nonzero)", mp.norm(Mq * Mb + Mb * Mq - 2 * s_ * mp.eye(4)))


# ------------------------------------------------------------------ N
def part_N():
    print("\nN. Neutral kernel in position space; boundary forms; the superspace flux; imaginary frequency")
    print(" The d=1 chain formula, in its convergent range")
    al = mp.mpf("0.3")
    # The tails fall as |z|^(-1.2), so map them to (0, 1) by z -> 1/u. That leaves u^(-0.8) at the origin,
    # which tanh-sinh resolves only to (precision floor)^0.2, so use high working precision.
    with mp.workdps(90):
        direct = (
            mp.quad(lambda z: z ** (-2 * al) * (1 - z) ** (-2 * al), [0, half, 1])
            + mp.quad(lambda u: u ** (4 * al - 2) * (1 - u) ** (-2 * al), [0, half, 1])
            + mp.quad(lambda z: z ** (-2 * al) * (1 + z) ** (-2 * al), [0, 1])
            + mp.quad(lambda u: u ** (4 * al - 2) * (1 + u) ** (-2 * al), [0, 1])
        )
    # the same four pieces in closed form
    record("the pieces are B(.4,.4) + 2 B(.2,.4)", mp.beta(1 - 2 * al, 1 - 2 * al) + 2 * mp.beta(4 * al - 1, 1 - 2 * al), direct, tol=mp.mpf(10) ** -8)
    cf = mp.sqrt(mp.pi) * mp.gamma(half - al) ** 2 * mp.gamma(2 * al - half) / (mp.gamma(al) ** 2 * mp.gamma(1 - 2 * al))
    record("int |z|^(-2a) |z-1|^(-2a) dz at a=0.3", direct, cf, tol=mp.mpf(10) ** -8)
    print(" The boson rail: int G_b G_b = -rho', with the same continuation")
    for q in (3, 5, 7):
        m = P(q)
        d, nu, mu = m.d, m.nu, m.mu
        al = (1 + d) / 2
        K = mp.sqrt(mp.pi) * mp.gamma(half - al) ** 2 * mp.gamma(2 * al - half) / (mp.gamma(al) ** 2 * mp.gamma(1 - 2 * al))
        C = 2 * mp.gamma(nu) * mp.gamma(2 * d - 1) / mp.gamma(nu + 2 * d - 1) - mp.beta(nu, nu)
        record(f"q={q}: d^2 K = 2 d rho_0 / b^2 = 4 d mu C", d**2 * K, 4 * d * mu * C)
    print(" First row on (f, -f'): K^11 f + K^12 (-f') against K^(0) f, smooth test functions, direct quadrature")
    for q in (3, 5):
        m = P(q)
        d, nu, mu = m.d, m.nu, m.mu
        b = mp.mpf(1)
        J = mp.tan(mp.pi * d / 2) / (2 * mp.pi)
        C = 2 * mp.gamma(nu) * mp.gamma(2 * d - 1) / mp.gamma(nu + 2 * d - 1) - mp.beta(nu, nu)
        Arail = lambda t: b**2 * C * abs(t) ** (2 * mu)
        rho = lambda t: 2 * mu * b**2 * C * mp.sign(t) * abs(t) ** (-2 * d)
        V = lambda y: J * (q - 1) * b ** (q - 2) * mp.sign(y) * abs(y) ** (-2 * mu)
        Gq2 = lambda y: b ** (q - 2) * mp.sign(y) * abs(y) ** (-2 * mu)  # G_psi^(q-2)
        # (q-2) G_b G_psi^(q-3), written out from G_b = d b / |y|^(1+d), G_psi^(q-3) = b^(q-3) |y|^(-d(q-3))
        GbG = lambda y: (q - 2) * (d * b / abs(y) ** (1 + d)) * b ** (q - 3) * abs(y) ** (-d * (q - 3))
        for name, f, df in (
            ("y^2 e^(-y^2)", lambda y: y**2 * mp.exp(-(y**2)), lambda y: (2 * y - 2 * y**3) * mp.exp(-(y**2))),
            ("y^3 e^(-y^2)", lambda y: y**3 * mp.exp(-(y**2)), lambda y: (3 * y**2 - 2 * y**4) * mp.exp(-(y**2))),
        ):
            for x in (mp.mpf("0.7"), mp.mpf("-1.3")):
                pts = sorted([-mp.inf, min(-x, mp.mpf(0)), max(-x, mp.mpf(0)), mp.inf])
                row1 = J * (q - 1) * mp.quad(lambda y: Arail(x + y) * (GbG(y) * f(y) - Gq2(y) * df(y)), pts)
                K0f = mp.quad(lambda y: rho(x + y) * V(y) * f(y), pts)
                record(f"q={q} f={name} x={x}: first row = K^(0) f", row1, K0f, tol=mp.mpf(10) ** -7)
    print(" Boundary forms: the conserved Wronskians, with coefficient 2 rho and the Dirac form")
    for q in (3, 5):
        m = P(q)
        d, nu, mu = m.d, m.nu, m.mu
        for rho_ in (nu, d):
            w = mp.mpf("1.37")
            A0, B0 = conn(rho_, 0, w)
            A1, B1 = conn(rho_, 1, w)
            for t0 in (mp.mpf("0.3"), mp.mpf("1.1")):
                f0 = lambda t: reg(rho_, 0, w, t)
                f1 = lambda t: reg(rho_, 1, w, t)
                W_ = f0(t0) * mp.diff(f1, t0) - f1(t0) * mp.diff(f0, t0)
                record(f"q={q} rho={mp.nstr(rho_, 3)} theta={t0}: f0 f1' - f1 f0' = -2 rho (A0 B1 - B0 A1)", W_, -2 * rho_ * (A0 * B1 - B0 * A1), tol=mp.mpf(10) ** -8)
        Om = mp.mpf("1.87")
        w = Om - half
        chi0, chi1 = dirac(m, Om, 0), dirac(m, Om, 1)
        a0, b0 = conn(nu, 0, w)[0], 2 * nu * conn(nu, 0, w)[1] / (w - mu)
        a1, b1 = conn(nu, 1, w)[0], 2 * nu * conn(nu, 1, w)[1] / (w - mu)
        for t0 in (mp.mpf("0.3"), mp.mpf("1.1")):
            val = (chi0(t0).T * S1 * G1 * chi1(t0))[0]
            record(f"q={q} theta={t0}: psi0^T s1 g1 psi1 = b0 a1 - a0 b1", val, b0 * a1 - a0 * b1, tol=mp.mpf(10) ** -8)
    print(" Superspace, with Grassmann variables: flux and superpotential on all components")
    import sympy as sp

    from check_supercasimir import odd

    th, tb, a, ah, bb, bh, da, dah = [odd(n) for n in ("t1", "tb1", "t2", "tb2", "t3", "tb3", "t4", "tb4")]
    A_, B_, Ab_, Bp_, dA_, dB_, c_ = sp.symbols("A B Ab Bp dA dB c")
    Pm = A_ + th * a + tb * ah + th * tb * B_
    dPm = dA_ + th * da + tb * dah + th * tb * dB_
    Pp = Ab_ + th * bb + tb * bh + th * tb * Bp_

    def top(x):
        return {k[2:]: v for k, v in x.terms.items() if k[:2] == ("t1", "tb1")}

    flux = top(Pp * dPm)
    # expected: A_b dB + B_+ dA - b dahat + bhat da
    want = {(): Ab_ * dB_ + Bp_ * dA_, ("t3", "tb4"): sp.Integer(-1), ("tb3", "t4"): sp.Integer(1)}
    ok = all(sp.simplify(flux.get(k, 0) - v) == 0 for k, v in want.items()) and set(flux) == set(want)
    RESULTS.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] top component of Phi_+ dPhi_-: {flux}")
    sq = top(Pm * Pm)
    ok = sp.simplify(sq.get((), 0) - 2 * A_ * B_) == 0 and sq.get(("t2", "tb2"), 0) == -2 and set(sq) == {(), ("t2", "tb2")}
    RESULTS.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] top component of Phi_-^2: {sq}")
    print(" Imaginary frequency: k_p(1/2 + iy) = s_p L (mu - iy) D_psi, k_(1-p)(-1/2 + iy) = s_p L (mu + iy) D_psi, D_psi real")
    for q in (3, 5):
        m = P(q)
        for p in (0, 1):
            sg = 1 if p == 0 else -1
            kp, kpb = (m.kS, m.kA) if p == 0 else (m.kA, m.kS)
            for y in (mp.mpf("0.3"), mp.mpf("1.7")):
                Dpsi = -1 / (2 * m.nu * resp(m.nu, p, I * y))
                record(f"q={q} p={p} y={y}: Im D_psi", mp.im(Dpsi), 0, tol=mp.mpf(10) ** -12)
                record(f"q={q} p={p} y={y}: k_p", kp(half + I * y), sg * m.L * (m.mu - I * y) * Dpsi, tol=mp.mpf(10) ** -10)
                record(f"q={q} p={p} y={y}: k_(1-p)", kpb(-half + I * y), sg * m.L * (m.mu + I * y) * Dpsi, tol=mp.mpf(10) ** -10)


# ------------------------------------------------------------------ O
def part_O():
    print("\nO. Fock space of one multiplet. Modes: b (boson at |w|, w = H - 1/2), f, fbar (Dirac, energy H), t (boson at H + 1/2)")
    Nb = 6

    def kron(*ms):
        out = ms[0]
        for x in ms[1:]:
            r, c = out.rows * x.rows, out.cols * x.cols
            new = mp.matrix(r, c)
            for i in range(out.rows):
                for j in range(out.cols):
                    if out[i, j] != 0:
                        for k in range(x.rows):
                            for l_ in range(x.cols):
                                new[i * x.rows + k, j * x.cols + l_] = out[i, j] * x[k, l_]
            out = new
        return out

    ab = mp.matrix(Nb, Nb)
    for n in range(1, Nb):
        ab[n - 1, n] = mp.sqrt(n)
    fm = mp.matrix([[0, 1], [0, 0]])
    Zf = mp.matrix([[1, 0], [0, -1]])
    Ib, If = mp.eye(Nb), mp.eye(2)
    # ordering: (N_b, N_f); and a separate block for (N_t, N_fbar) of the same form
    aB = kron(ab, If)
    aF = kron(Ib, fm)
    for H in (mp.mpf("0.3"), mp.mpf("0.4"), mp.mpf("0.45"), mp.mpf("0.7"), mp.mpf("1.3")):
        w = H - half
        NB, NF = aB.T * aB, aF.T * aF
        # the charge lowering L_0 by 1/2: for w > 0 it trades f for b; for w < 0 the boson image mode
        # e^(-i w t) = e^(+i|w| t) carries a creation operator in Phi, so the charge is a_b a_f
        Qbar = mp.sqrt(2 * abs(w)) * (aB * aF if w < 0 else aB.T * aF)
        Q = Qbar.T
        An = Q * Qbar + Qbar * Q
        Evac = (abs(w) + H + half) / 2 - H
        # L_0 - J_0/2 on this pair; the (fbar, t) pair is positive and contributes (H + 1/2)(N_t + N_fbar) separately
        Htw = Evac * mp.eye(2 * Nb) + abs(w) * NB + (H - half) * NF
        low = [i for i in range(2 * Nb) if i // 2 < Nb - 2]
        res = max(abs(An[i, j] - 2 * Htw[i, j]) for i in low for j in low)
        record(f"H={H}: {{Q, Qbar}} = 2 (L_0 - J_0/2) on low states, with E_vac from the zero-point sum", res, 0, tol=mp.mpf(10) ** -15)
        vac = 0  # index of |N_b=0, N_f=0>
        filled = 1  # index of |N_b=0, N_f=1>
        normQ = sum(abs(Q[i, vac]) ** 2 + abs(Qbar[i, vac]) ** 2 for i in range(2 * Nb))
        normF = sum(abs(Q[i, filled]) ** 2 + abs(Qbar[i, filled]) ** 2 for i in range(2 * Nb))
        record(f"H={H}: |Q vac|^2 + |Qbar vac|^2 = 2 E_vac", normQ, 2 * Evac, tol=mp.mpf(10) ** -15)
        if w < 0:
            record(f"H={H}: the state with f filled is invariant", normF, 0, tol=mp.mpf(10) ** -15)
            info(f"   that state: L_0 = E_vac + H, J_0 = 1, L_0 - J_0/2", f"{mp.nstr(Evac + H, 6)}, 1, {mp.nstr(Evac + H - half, 6)}")
    _ = Zf


# ------------------------------------------------------------------ P
def part_P():
    print("\nP. Bulycheva (4.99): chi^h read as chi^h 2F1(h,h;2h+1;chi), summed over h = 2 Delta + n")
    print(" Printed factor: sin(pi h) - sin(2 pi D)/cos(pi h). With (4.39) corrected: sin(2 pi D) - sin(pi h)/cos(pi h).")
    for q in (3, 5):
        Dl = mp.mpf(1) / (2 * q)
        d = 2 * Dl
        printed = lambda h: mp.sin(mp.pi * h) - mp.sin(mp.pi * d) / mp.cos(mp.pi * h)
        corrected = lambda h: mp.sin(mp.pi * d) - mp.sin(mp.pi * h) / mp.cos(mp.pi * h)

        def coef(n, fac):
            h = d + n
            # -Res_h Gamma(2 Delta - h) at h = 2 Delta + n is (-1)^n / n!
            return (
                (-1) ** n / mp.factorial(n)
                * (1 / (8 * mp.pi * mp.tan(mp.pi * h)))
                * (16 * Dl / mp.tan(mp.pi * Dl))
                * mp.gamma(h) ** 2 * mp.gamma(d + h) * mp.gamma(-d) / (mp.gamma(2 * h) * mp.gamma(d))
                * fac(h)
            )

        record(f"q={q}: leading coefficient of (4.99) as printed", coef(0, printed), 1, tol=mp.mpf(10) ** -12)
        info(f"   printed / corrected at n = 0, 1, 2, 3", ", ".join(mp.nstr(printed(d + n) / corrected(d + n), 5) for n in range(4)))
        # The coefficient of phi_h in k(h) xi_h + k(-h) xi_(-h), with k(h) = K (s - sin pi h), K even in h
        h = mp.mpf("0.37")
        sec = 1 / mp.cos(mp.pi * h)
        S, s_ = mp.sin(mp.pi * h), mp.sin(mp.pi * d)
        record(f"q={q}: corrected (4.39) gives the corrected factor", (s_ - S) * (1 + sec) + (s_ + S) * (1 - sec), 2 * corrected(h))
        record(f"q={q}: printed (4.39) gives the printed factor", (s_ - S) * (1 + sec) - (s_ + S) * (1 - sec), -2 * printed(h))
        for chi in (mp.mpf("0.2"), mp.mpf("0.5"), mp.mpf("0.8")):
            total = mp.nsum(lambda n: coef(int(n), corrected) * chi ** (d + int(n)) * mp.hyp2f1(d + int(n), d + int(n), 2 * (d + int(n)) + 1, chi), [0, mp.inf])
            record(f"q={q} chi={chi}: sum with the corrected factor and the block = chi^(2 Delta)", total, chi**d, tol=mp.mpf(10) ** -8)
            bad = mp.nsum(lambda n: coef(int(n), printed) * chi ** (d + int(n)) * mp.hyp2f1(d + int(n), d + int(n), 2 * (d + int(n)) + 1, chi), [0, mp.inf])
            bare = mp.nsum(lambda n: coef(int(n), printed) * chi ** (d + int(n)), [0, mp.inf])
            info(f"   printed factor, with the block / with bare chi^h, for comparison", f"{mp.nstr(bad, 8)} / {mp.nstr(bare, 8)}")


PARTS = {"A": part_A, "B": part_B, "C": part_C, "D": part_D, "E": part_E, "F": part_F, "G": part_G, "H": part_H, "I": part_I, "J": part_J, "K": part_K, "M": part_M, "e": part_E2, "N": part_N, "O": part_O, "P": part_P}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dps", type=int, default=20)
    parser.add_argument("--parts", default="ABCDEeFGHIJKMNOP")
    args = parser.parse_args()
    mp.mp.dps = args.dps
    for key in args.parts:
        PARTS[key]()
    failed = RESULTS.count(False)
    print(f"\n{len(RESULTS) - failed}/{len(RESULTS)} checks passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
