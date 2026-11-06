#!/usr/bin/env python3
"""Numerical checks of the Robin-ladder identities in Thong (2026) and the N=2 draft.

Tests closed-form statements from
  [T]  K. H. Thong, arXiv:2610.31989v1 (papers/thong_2026.pdf), equations (x.y);
  [N2] research/pdf/n2_susy_generalization_note.pdf, equations (n).

The script tests internal consistency of the printed formulas. It does NOT check the
N=2 component kernels against Peng-Spradlin-Volovich, which are not in papers/.

Usage:
    python scripts/check_robin_identities.py [--dps 30]

Exits non-zero if any check fails. Memory use is negligible. No random numbers are used; all
test points are fixed in the code. Results of a recorded run: results/data/robin_identity_checks/.
"""

import argparse
import sys

import mpmath as mp

RESULTS = []


def record(label, value, reference, tol=None):
    """Record |value - reference| against a tolerance."""
    tol = mp.mpf(10) ** (-mp.mp.dps // 2) if tol is None else tol
    residual = abs(mp.mpc(value) - mp.mpc(reference))
    passed = residual <= tol
    RESULTS.append(passed)
    print(f"[{'PASS' if passed else 'FAIL'}] {label}: residual={mp.nstr(residual, 3)}")


def bisect_root(f, lo, hi):
    return mp.findroot(f, (mp.mpf(lo), mp.mpf(hi)), solver="anderson")


# ---------------------------------------------------------------------------
# Common response function: [T] (B.2), [N2] (6). At d=1, J labels strip parity.
# ---------------------------------------------------------------------------
def response(rho, J, omega, d=1):
    plus, minus = mp.mpf(d) / 2 + rho, mp.mpf(d) / 2 - rho
    value = mp.gamma(-rho) / mp.gamma(rho)
    for sign in (1, -1):
        value *= mp.gamma((J + plus + sign * omega) / 2)
        value /= mp.gamma((J + minus + sign * omega) / 2)
    return value


# ---------------------------------------------------------------------------
# Thong (2026)
# ---------------------------------------------------------------------------
def kc_syk(q, h):
    """Conformal SYK kernel eigenvalue, Maldacena-Stanford eq. (3.35).

    Transcribed from memory: that paper is not in papers/ yet. kc(2)=1 and Thong's relation
    kc(h) R_-(h-1/2) = 4^-nu (q-1) both hold with it, but compare it to the source once available.
    """
    d = mp.mpf(1) / q
    prefactor = (
        -(q - 1)
        * mp.gamma(1.5 - d)
        * mp.gamma(1 - d)
        / (mp.gamma(0.5 + d) * mp.gamma(d))
    )
    return (
        prefactor
        * mp.gamma(d + h / 2)
        * mp.gamma(0.5 + d - h / 2)
        / (mp.gamma(1.5 - d - h / 2) * mp.gamma(1 - d + h / 2))
    )


def check_thong():
    print("\n== Thong (2026) ==")
    for q in (4, 6, 10):
        nu = 1 - mp.mpf(2) / q
        kappa = mp.power(4, -nu) * (q - 1)
        for h in (mp.mpf("2.7"), mp.mpc("1.3", "0.8")):
            record(
                f"(C.1 text) kc(h) R_-(h-1/2) = 4^-nu (q-1), q={q}, h={h}",
                kc_syk(q, h) * response(nu, 1, h - 0.5),
                kappa,
            )
        record(f"kc(2)=1, q={q}", kc_syk(q, mp.mpf(2)), 1)
        record(
            f"Robin root at h=2 is omega=3/2, q={q}", response(nu, 1, mp.mpf(1.5)), kappa
        )
        record(
            f"(E.3) DGJS coupling = 4^-nu (q-1), q={q}",
            (q - 1)
            * mp.gamma(1.5 - mp.mpf(1) / q)
            * mp.gamma(1 - mp.mpf(1) / q)
            * mp.gamma(mp.mpf(2) / q)
            / (
                mp.gamma(0.5 + mp.mpf(1) / q)
                * mp.gamma(2 - mp.mpf(2) / q)
                * mp.gamma(mp.mpf(1) / q)
            ),
            kappa,
        )
        record(
            f"(C.3, A.6 text) c_boundary = -2 nu 4^-nu (q-1), q={q}",
            cboundary_syk(q),
            -2 * nu * mp.power(4, -nu) * (q - 1),
        )

    # Bilinear spectrum at q=4: the Robin roots against direct roots of kc=1.
    nu = mp.mpf(1) / 2
    kappa = mp.mpf(3) / 2
    print("     q=4 spectrum from R_-(omega)=kappa, h = omega + 1/2:")
    for lo, hi in ((3.51, 3.99), (5.51, 5.99), (7.51, 7.99)):
        h_robin = bisect_root(lambda h: response(nu, 1, h - 0.5) - kappa, lo, hi)
        h_kernel = bisect_root(lambda h: kc_syk(4, h) - 1, lo, hi)
        print(f"       h = {mp.nstr(h_robin, 12)}")
        record(f"q=4 root near {lo}", h_robin, h_kernel)

    # Neumann identification (3.4)-(3.5): cylinder kernel against -1/(2 nu R_J).
    for nu_test in (mp.mpf("0.3"), mp.mpf("0.8")):
        a = mp.mpf(1) / 2 - nu_test
        for p in (mp.mpf("0.7"), mp.mpf("2.1")):
            prefactor = mp.gamma(a) / (2 * mp.sqrt(mp.pi) * mp.gamma(1 - nu_test))
            # 2 cosh t -+ 2 written as 4 sinh^2(t/2), 4 cosh^2(t/2) to avoid cancellation.
            integrand = lambda t: (
                (4 * mp.sinh(t / 2) ** 2) ** (-a) - (4 * mp.cosh(t / 2) ** 2) ** (-a)
            ) * mp.cos(p * t) if t != 0 else mp.mpf(0)
            integral = 2 * mp.quad(integrand, [0, 1, 4, 12, 40])
            record(
                f"(3.4) d=1 odd Neumann kernel, nu={nu_test}, p={p}",
                prefactor * integral,
                -1 / (2 * nu_test * response(nu_test, 1, 1j * p)),
                tol=mp.mpf("1e-8"),
            )
    # d=2, J=0: the AR parameters nu=2/3, a=1/3.
    nu_ar, a_ar = mp.mpf(2) / 3, mp.mpf(1) / 3
    for p in (mp.mpf("0.7"), mp.mpf("1.6")):
        # Angular integral in closed form; the tail decays as exp(-2 a t) = exp(-2t/3).
        angular = lambda t: (
            2 * mp.pi * (2 * mp.cosh(t)) ** (-a_ar)
            * mp.hyp2f1(a_ar / 2, (a_ar + 1) / 2, 1, mp.sech(t) ** 2)
        )
        breakpoints = [0, 1] + [1 + 4 * k for k in range(1, 26)]
        integral = 2 * mp.quad(lambda t: angular(t) * mp.cos(p * t), breakpoints)
        record(
            f"(3.4) d=2 J=0 Neumann kernel (AR), p={p}",
            mp.gamma(a_ar) / (2 * mp.pi * mp.gamma(1 - nu_ar)) * integral,
            -1 / (2 * nu_ar * response(nu_ar, 0, 1j * p, d=2)),
            tol=mp.mpf("1e-5"),
        )

    # (C.6) centre integral of a block.
    for h in (mp.mpf("2.4"), mp.mpf("3.1")):
        eta = mp.mpf("1.7")
        block = lambda chi: chi**h * mp.hyp2f1(h, h, 2 * h, chi)
        lhs = 2 * mp.quad(lambda u: block(2 / (1 + eta + 2 * u**2)), [0, 1, 5, mp.inf])
        omega = h - 0.5
        norm = mp.gamma(h) ** 2 / (2 * mp.gamma(2 * h))
        record(
            f"(C.6) J_h(eta), h={h}",
            lhs,
            mp.pi / (2 * omega * norm) * mp.exp(-omega * mp.acosh(eta)),
            tol=mp.mpf("1e-10"),
        )

    # AR threshold (4.7) and normalisation of appendix D.
    gc = -8 * mp.pi / 3 * response(nu_ar, 0, mp.mpf(0), d=2)
    print(f"     AR g_c = {mp.nstr(gc, 12)} (paper: 1.02236136)")
    record("(4.7) AR g_c", gc, mp.mpf("1.02236136"), tol=mp.mpf("1e-8"))
    b = mp.gamma(mp.mpf(2) / 3) / (mp.power(2, mp.mpf(2) / 3) * mp.pi * mp.gamma(mp.mpf(1) / 3))
    h_phi = b**2 * mp.pi * mp.gamma(mp.mpf(1) / 3) ** 3 / mp.gamma(mp.mpf(2) / 3) ** 3
    record(
        "(A.5),(A.6),(D.1) AR c_boundary = 1/(2 pi)",
        4 * mp.pi * (b / 2) * h_phi,
        1 / (2 * mp.pi),
    )

    for q in (4, 6, 10):
        nu = 1 - mp.mpf(2) / q
        assert mp.power(4, -nu) * (q - 1) > response(nu, 1, mp.mpf(0))
    print("     SYK kappa > R_-(0) for q=4,6,10: stability criterion of app. B holds")


def cboundary_syk(q):
    """c_boundary from the chain integral (C.3), (A.6), with epsilon=-1 from (C.1).

    The Gamma(a) in the denominator of (A.6) is singular at q=4 (a=0). The bracket in
    (C.3) vanishes there too, so the q=4 value is taken as a limit in Delta.
    """
    if q == 4:
        return mp.limit(lambda eps: cboundary_syk(4 + eps), 0, direction=1)
    d = mp.mpf(1) / q
    nu = 1 - 2 * d
    a = 2 * d - mp.mpf(1) / 2
    j2bq = (mp.mpf(1) / 2 - d) * mp.tan(mp.pi * d) / mp.pi
    bracket = 2 * mp.beta(nu, 1 - 2 * nu) - mp.beta(nu, nu)
    return -2 * mp.sqrt(mp.pi) * (q - 1) * j2bq * bracket * mp.gamma(1 - nu) / mp.gamma(a)


# ---------------------------------------------------------------------------
# N=2 draft
# ---------------------------------------------------------------------------
class N2:
    def __init__(self, q):
        self.q = q
        self.delta = mp.mpf(1) / q
        self.nu = 1 - self.delta
        self.mu = mp.mpf(1) / 2 - self.delta
        self.L = 2 * (q - 1) * mp.power(4, -self.nu)
        self.C = self.nu * mp.gamma(self.nu) / (mp.pi * self.delta * mp.gamma(self.delta))
        self.norm = self.delta / (2 * self.nu * mp.tan(mp.pi * self.delta / 2))
        self.D = 2 * self.norm / mp.pi

    def kA(self, H):  # (15)
        d = self.delta
        return (
            -self.C
            * mp.gamma(d - H)
            * mp.gamma(d + H)
            * (mp.sin(mp.pi * d) - mp.sin(mp.pi * H))
        )

    def kS(self, H):  # (15)
        d = self.delta
        return (
            -self.C
            * mp.gamma(d - H)
            * mp.gamma(d + H)
            * (mp.sin(mp.pi * d) + mp.sin(mp.pi * H))
        )

    def sigma(self, p):
        return 1 if p == 0 else -1

    def Rnu(self, p, omega):
        return response(self.nu, p, omega)

    def Rdelta(self, p, omega):
        return response(self.delta, 1 - p, omega)

    def Dmat(self, p, omega):  # (7)
        return mp.matrix(
            [
                [-1 / (2 * self.nu * self.Rnu(p, omega)), 0],
                [0, self.Rdelta(p, omega) / (2 * self.delta)],
            ]
        )

    def Mmat(self, p):  # (8)
        factor = self.sigma(p) * self.L
        return mp.matrix(
            [[factor * 2 * self.mu, factor * 2 * self.delta], [factor * 2 * self.delta, 0]]
        )

    def kbar(self, p, omega):  # (9)
        return self.Dmat(p, omega) * self.Mmat(p)

    def phi(self, H, z):  # (24)
        return (
            mp.gamma(H) ** 2 / mp.gamma(2 * H) * z**H * mp.hyp2f1(H, H, 2 * H + 1, z)
        )

    def fbos(self, H, z):
        return mp.gamma(H) ** 2 / mp.gamma(2 * H) * z**H * mp.hyp2f1(H, H, 2 * H, z)

    def contour(self, weight, z):
        """Integrate upward through Re H = delta/2, tilted right to gain decay.

        Valid when all poles of `weight` are on the real axis to the right of the
        crossing point, which is the case at order g^0 and g^1.
        """
        c = self.delta / 2
        phase = mp.expjpi(mp.mpf(1) / 3)
        f = lambda t: weight(c + t * phase) * self.phi(c + t * phase, z) * phase
        return mp.im(mp.quad(f, [0, 2, 8, 30, 80])) / mp.pi

    def seed_weight(self, H):  # order g^0 of (26)
        return self.norm * (
            mp.cot(mp.pi * H / 2) * self.kA(H) - mp.tan(mp.pi * H / 2) * self.kS(H)
        )

    def first_rung_weight(self, H):  # order g^1 of (26)
        return self.norm * (
            mp.cot(mp.pi * H / 2) * self.kA(H) ** 2
            - mp.tan(mp.pi * H / 2) * self.kS(H) ** 2
        )

    def first_rung_direct(self, z):  # (34)
        d, nu = self.delta, self.nu
        P = nu * (1 - 2 * d) / (2 * mp.pi * d) * mp.tan(mp.pi * d / 2)
        Cconv = 2 * mp.beta(nu, 2 * d - 1) - mp.beta(nu, nu)

        def integrand(v):
            # Rounding can land a quadrature node on an integrable endpoint singularity,
            # where the node's weight is negligible; drop such nodes.
            if v == 0 or v == z or v == 1:
                return mp.mpf(0)
            return (
                mp.sign(v)
                * mp.sign(z - v)
                * mp.sign(v - 1)
                * abs(v) ** (-nu)
                * abs(z - v) ** (-d)
                * abs(v - 1) ** (-d)
            )

        # v = +-u^q removes the slow |v|^(-nu) endpoint: dv |v|^(-nu) = q du at v -> 0.
        near_zero = lambda sign, upper: mp.quad(
            lambda u: self.q * u ** (self.q - 1) * integrand(sign * u**self.q),
            [0, upper ** (mp.mpf(1) / self.q)],
        )
        integral = (
            mp.quad(integrand, [-mp.inf, -1])
            + near_zero(-1, mp.mpf(1))
            + near_zero(1, z)
            + mp.quad(integrand, [z, 1, 2, mp.inf])
        )
        return P * Cconv * z**d * integral


def check_n2():
    print("\n== N=2 draft ==")
    for q in (3, 5, 7):
        m = N2(q)
        d, nu, mu = m.delta, m.nu, m.mu
        record(f"(5) m_psi^2 = nu^2-1/4, q={q}", mu * (mu + 1), nu**2 - 0.25)
        record(f"(5) m_b^2 = delta^2-1/4, q={q}", mu * (mu - 1), d**2 - 0.25)
        record(f"(15) kA(1)=1, q={q}", m.kA(mp.mpf(1)), 1)
        record(f"(15) kS(1)=1, q={q}", m.kS(mp.mpf(1)), 1)
        record(f"(15) kA(0)=-(q-1), q={q}", m.kA(mp.mpf(0)), -(q - 1))

        for p in (0, 1):
            for h in (mp.mpf("2.3"), mp.mpc("1.7", "0.6")):
                omega = h - 0.5
                record(
                    f"(13) R_nu R_delta, q={q}, p={p}, h={h}",
                    m.Rnu(p, omega) * m.Rdelta(p, omega),
                    (mu**2 - omega**2) / (4 * d * nu),
                )
                kb = m.kbar(p, omega)
                first, second = (m.kA(h), m.kS(h - 1)) if p == 1 else (m.kS(h), m.kA(h - 1))
                record(
                    f"(9),(16) trace of k-hat, q={q}, p={p}, h={h}",
                    kb[0, 0] + kb[1, 1],
                    first + second,
                )
                record(
                    f"(9),(16) det of k-hat, q={q}, p={p}, h={h}",
                    kb[0, 0] * kb[1, 1] - kb[0, 1] * kb[1, 0],
                    first * second,
                )
                Omega = h
                RF = 2 * nu * m.Rnu(p, Omega - 0.5) / (Omega - nu)
                RF_top = 2 * nu * m.Rnu(1 - p, Omega + 0.5) / (Omega + 0.5 + mu)
                branch = m.kA(Omega) if p == 1 else m.kS(Omega)
                record(f"(22) k = sigma L/R_F, q={q}, p={p}, H={h}", branch, m.sigma(p) * m.L / RF)
                record(f"(sec.4) R_F top-scalar form, q={q}, p={p}, H={h}", RF, RF_top)

        z = mp.mpf("0.37")
        for H in (mp.mpf("0.9"), mp.mpf("2.2")):
            record(
                f"(24) phi = F_H - F_(H+1), q={q}, H={H}",
                m.phi(H, z),
                m.fbos(H, z) - m.fbos(H + 1, z),
            )
        for z in (mp.mpf("0.2"), mp.mpf("0.6")):
            record(
                f"(29) seed identity, q={q}, z={z}",
                mp.nsum(
                    lambda n: (-1) ** int(n)
                    * mp.gamma(2 * d + n)
                    / (mp.gamma(d) ** 2 * mp.factorial(n))
                    * m.phi(d + n, z),
                    [0, mp.inf],
                ),
                z**d,
                tol=mp.mpf("1e-10"),
            )

        # (12): fermion rail against the chain integral, and its ratio to the cylinder
        # coefficient. Both are divided by b_psi^2, using J b^q = tan(pi delta/2)/(2 pi).
        rail = (
            -2 * mp.gamma(nu) ** 2 * mp.cos(mp.pi * d / 2) ** 2
            / (mp.gamma(2 * nu) * mp.cos(mp.pi * d))
        )
        record(
            f"(12) r_psi = chain integral, q={q}",
            rail,
            2 * mp.beta(nu, 2 * d - 1) - mp.beta(nu, nu),
        )
        c_psi = mp.gamma(d - 0.5) / (2 * mp.sqrt(mp.pi) * mp.gamma(d))
        record(
            f"(12) r_psi/c_psi = l_psi, q={q}",
            rail / c_psi,
            2 * mp.power(4, -nu) * 2 * mp.pi / mp.tan(mp.pi * d / 2),
        )

        # The ordinary Abel map behind (25), inverse to Thong's (C.6): A exp(-(H-1/2)t) = F_H/2.
        for H, z in ((mp.mpf("0.9"), mp.mpf("0.4")), (mp.mpf("2.2"), mp.mpf("0.7"))):
            eta = mp.acosh(2 / z - 1)
            abel = lambda s: (
                (1 / mp.sqrt(2)) * 2 * mp.sqrt(2) * mp.exp(-(H - 0.5) * (eta + s**2))
                / mp.sqrt(2 * mp.sinh(eta + s**2 / 2) * _sinhc(s**2 / 2))
            )
            record(
                f"(25) A exp(-(H-1/2)t) = F_H/2, q-independent, H={H}, z={z}",
                mp.quad(abel, [0, 1, 3, 12]),
                m.fbos(H, z) / 2,
                tol=mp.mpf("1e-10"),
            )

        b1 = mp.diff(m.kS, 1)
        record(f"(31) b1, q={q}", b1, (1 - 2 * d) / (d * nu) - mp.pi / mp.sin(mp.pi * d), tol=mp.mpf("1e-10"))
        record(
            f"(31) b2, q={q}",
            mp.diff(m.kS, 1, 2),
            b1**2 + mp.polygamma(1, d) - mp.polygamma(1, nu) + nu**-2 - d**-2,
            tol=mp.mpf("1e-8"),
        )
        record(
            f"(sec.6) kA'(0), q={q}",
            mp.diff(m.kA, mp.mpf("1e-30")),
            mp.pi * nu / (d * mp.sin(mp.pi * d)),
            tol=mp.mpf("1e-10"),
        )

    m = N2(3)
    H_S = bisect_root(lambda H: m.kS(H) - 1, 2.5, 3.2)
    H_A = bisect_root(lambda H: m.kA(H) - 1, 2.6, 3.3)
    record("(sec.3) H_S, q=3", H_S, mp.mpf("2.8211441710"), tol=mp.mpf("1e-9"))
    record("(sec.3) H_A, q=3", H_A, mp.mpf("3.0659046429"), tol=mp.mpf("1e-9"))
    b1 = mp.diff(m.kS, 1)
    b2 = mp.diff(m.kS, 1, 2)
    record("(33) coefficient 1, q=3", m.D / b1, mp.mpf("-0.129565995703"), tol=mp.mpf("1e-11"))
    record(
        "(33) coefficient 2, q=3",
        m.D * (1 - b2 / (2 * b1**2)),
        mp.mpf("0.129254113198"),
        tol=mp.mpf("1e-10"),
    )

    # The Abel transform (25), explicit form on the page.
    for H, z in ((mp.mpf("0.8"), mp.mpf("0.4")), (mp.mpf("2.3"), mp.mpf("0.6"))):
        eta = mp.acosh(2 / z - 1)
        # t = eta + s^2 and cosh t - cosh eta = 2 sinh(eta + s^2/2) sinh(s^2/2) make the
        # inverse-square-root endpoint regular: 2s/sqrt(sinh(s^2/2)) -> 2 sqrt(2).
        integrand = lambda s: (
            2 * mp.sqrt(2) * mp.sinh((eta + s**2) / 2)
            * mp.exp(-H * (eta + s**2))
            / mp.sqrt(2 * mp.sinh(eta + s**2 / 2) * _sinhc(s**2 / 2))
        )
        record(
            f"(25) T2 exp(-Ht) = phi_H, H={H}, z={z}",
            2 * mp.sqrt(2) * mp.quad(integrand, [0, 1, 3, 10]),
            m.phi(H, z),
            tol=mp.mpf("1e-10"),
        )

    # Spectral reconstruction (27) at low order against the seed and direct convolution (34).
    table = {
        (3, "0.2"): mp.mpf("-0.533994288253294417"),
        (3, "0.6"): mp.mpf("0.130574305406873014"),
        (5, "0.4"): mp.mpf("-0.371415562926934307"),
    }
    for (q, z_text), printed in table.items():
        m = N2(q)
        z = mp.mpf(z_text)
        direct = m.first_rung_direct(z)
        record(f"(34) direct convolution vs table, q={q}, z={z_text}", direct, printed, tol=mp.mpf("1e-10"))
        spectral = m.contour(m.first_rung_weight, z) + m.D * mp.nsum(
            lambda n: m.kA(mp.mpf(int(n))) ** 2 * m.phi(mp.mpf(int(n)), z), [1, mp.inf]
        )
        record(f"d/dg of (27) at g=0 vs table, q={q}, z={z_text}", spectral, printed, tol=mp.mpf("1e-8"))
        seed = m.contour(m.seed_weight, z) + m.D * mp.nsum(
            lambda n: m.kA(mp.mpf(int(n))) * m.phi(mp.mpf(int(n)), z), [1, mp.inf]
        )
        record(f"(27) at g=0 reproduces z^delta, q={q}, z={z_text}", seed, z**m.delta, tol=mp.mpf("1e-8"))


def _sinhc(x):
    return mp.sinh(x) / x if x != 0 else mp.mpf(1)


def check_seed_mode_sum():
    """Test the mode content of the reduced seed kernel (3.4), including a + J < 0.

    Claim: Gamma(a)/(2 sqrt(pi) Gamma(1-nu)) * (D_s + sigma D_o) = sum_n w_n exp(-omega_n |tau|)
    with omega_n = a + J + 2n and w_n = 1/(2 nu R_J'(omega_n)), where the n=0 term keeps its
    continued exponent. When a + J < 0 this term grows. The stationary covariance, with the
    same Fourier response, has w(|a+J|) exp(-|a+J| |tau|), and the difference is
    -2 w(|a+J|) cosh((a+J) tau), which solves the homogeneous equation.
    """
    print("\n== Reduced seed: mode sum and the a + J < 0 sector ==")
    for q in (3, 5):
        m = N2(q)
        nu, a = m.nu, -m.mu
        prefactor = mp.gamma(a) / (2 * mp.sqrt(mp.pi) * mp.gamma(1 - nu))
        for J, sigma in ((0, 1), (1, -1)):
            weight = lambda w: 1 / (2 * nu * mp.diff(lambda x: response(nu, J, x), w))
            for tau in (mp.mpf("0.9"), mp.mpf("2.1")):
                kernel = prefactor * (
                    (4 * mp.sinh(tau / 2) ** 2) ** (-a)
                    + sigma * (4 * mp.cosh(tau / 2) ** 2) ** (-a)
                )
                modes = mp.nsum(
                    lambda n: weight(a + J + 2 * int(n)) * mp.exp(-(a + J + 2 * int(n)) * tau),
                    [0, mp.inf],
                )
                record(
                    f"(3.4) seed mode sum, q={q}, J={J}, a+J={mp.nstr(a + J, 4)}, tau={tau}",
                    kernel,
                    modes,
                    tol=mp.mpf("1e-10"),
                )
            if a + J < 0:
                w_stationary = weight(-(a + J))
                record(
                    f"w(-|a+J|) = -w(|a+J|), q={q}, J={J}",
                    weight(a + J),
                    -w_stationary,
                    tol=mp.mpf("1e-10"),
                )
                print(
                    f"     q={q} J={J}: continued seed has {mp.nstr(-w_stationary, 6)} "
                    f"exp(+{mp.nstr(-(a + J), 4)}|tau|); stationary covariance has "
                    f"{mp.nstr(w_stationary, 6)} exp(-{mp.nstr(-(a + J), 4)}|tau|)"
                )


def _positive_roots(m, p, g, omega_max):
    """Positive roots of det(I - g k-hat) from the scalar branches, by sign changes."""
    branches = (
        (lambda w: m.kA(w + 0.5), lambda w: m.kS(w - 0.5))
        if p == 1
        else (lambda w: m.kS(w + 0.5), lambda w: m.kA(w - 0.5))
    )
    names = ("h-branch", "(h-1)-branch")
    roots = []
    for branch, name in zip(branches, names):
        grid = [mp.mpf(i) * mp.mpf("0.0125") + mp.mpf("0.00137") for i in range(int(omega_max / 0.0125))]
        f = lambda w: mp.re(branch(w)) - 1 / mp.mpf(g)
        values = [(w, f(w)) for w in grid]
        for (w0, f0), (w1, f1) in zip(values, values[1:]):
            # A sign change with small values is a root; a large one is a pole.
            if f0 * f1 < 0 and abs(f0) < 50 and abs(f1) < 50:
                try:
                    roots.append((bisect_root(f, w0, w1), name))
                except (ValueError, ZeroDivisionError):
                    pass
    return sorted(roots, key=lambda item: item[0])


def diagnose_n2_modes():
    """Not a check of the draft: maps auxiliary modes against Thong's validity conditions."""
    print("\n== N=2 auxiliary-mode diagnostic ==")
    for q in (3, 5, 7):
        m = N2(q)
        gstar = 1 / m.kS(mp.mpf(1) / 2)
        print(
            f"q={q}: a_psi+J = {mp.nstr(-m.mu, 5)} (p=0), {mp.nstr(1 - m.mu, 5)} (p=1); "
            f"kS(1/2)={mp.nstr(m.kS(mp.mpf(1) / 2), 6)}, g*={mp.nstr(gstar, 8)}; "
            f"kA(1/2)={mp.nstr(m.kA(mp.mpf(1) / 2), 6)} (negative: no p=1 crossing at g>0)"
        )
        for g in (mp.mpf(1) / 4, gstar, mp.mpf(3) / 4, mp.mpf(1)):
            for p in (0, 1):
                roots = _positive_roots(m, p, g, 4.0)
                text = ", ".join(f"{mp.nstr(w, 6)} [{n}]" for w, n in roots[:4])
                print(f"   g={mp.nstr(g, 5)} p={p}: positive omega roots: {text}")
        for p in (0, 1):
            for g in (mp.mpf(1), mp.mpf("0.9"), mp.mpf("0.5"), mp.mpf("0.25")):
                for w, name in _positive_roots(m, p, g, 4.0)[:4]:
                    print(
                        f"   weight q={q} p={p} g={mp.nstr(g, 3)} omega={mp.nstr(w, 6)} "
                        f"[{name}]: eigenvalues of -Res G = {_residue_eigs(m, p, g, w)}"
                    )


def _residue_eigs(m, p, g, omega0):
    """Eigenvalues of minus the residue of G_{g,p} = (D^{-1} - g M)^{-1} at omega0."""
    M = m.Mmat(p)

    def G(w):
        D = m.Dmat(p, w)
        return mp.inverse(mp.inverse(D) - g * M)

    radius = mp.mpf("1e-3")
    residue = mp.matrix(2, 2)
    for k in range(32):
        phase = mp.expjpi(mp.mpf(2 * k + 1) / 32)
        value = G(omega0 + radius * phase)
        residue += value * radius * phase / 32
    residue = -residue.apply(mp.re)
    symmetrised = (residue + residue.T) / 2
    eigs = mp.eigsy(symmetrised, eigvals_only=True)
    return [mp.nstr(e, 5) for e in eigs]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dps", type=int, default=30, help="mpmath working decimal digits")
    args = parser.parse_args()
    mp.mp.dps = args.dps
    check_thong()
    check_n2()
    check_seed_mode_sum()
    diagnose_n2_modes()
    failures = RESULTS.count(False)
    print(f"\n{len(RESULTS) - failures}/{len(RESULTS)} checks passed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
