#!/usr/bin/env python3
"""Compare the N=2 draft's physical inputs with the formulas printed in its sources.

scripts/check_robin_identities.py shows the draft is internally consistent. This script
asks the different question of whether its inputs agree with the literature. Every
source formula below is transcribed from the PDF in papers/, with its equation number:

    [PSV]  Peng, Spradlin, Volovich, arXiv:1706.06078v2
    [B18]  Bulycheva, arXiv:1801.09006v2
    [FGMS] Fu, Gaiotto, Maldacena, Sachdev, arXiv:1610.08917
    [MS]   Maldacena, Stanford, arXiv:1604.07818v1 (v1 numbering, see the note on it)
    [T]    Thong, arXiv:2610.31989v1
    [B17]  Bulycheva, arXiv:1706.07411v1
    [GKST] Gu, Kitaev, Sachdev, Tarnopolsky, arXiv:1910.14099v3

Parts:
    1. Component kernels: [PSV] and [FGMS] against [N2](9), (15), (16), (22).
    2. The printed spectra of [PSV] against roots of the draft's kernels.
    3. [B18] kernel, including the sign the draft says is misprinted.
    4. [B18] Eq. (4.39), the other claimed misprint, and the seed normalization.
    5. [MS] kernel against kc_syk in check_robin_identities.py.
    6. Component structure of the superfield eigenfunctions.
    7. [PSV] Eqs. (4.10)-(4.26) plus [T](C.8) against [N2](12): the reduced four-point
       function at finite g, including the H=1 modes.
    8. The centre integral of the blocks, in the odd sector as a control on [T](C.8), and
       in the even sector, which [T] does not treat.
    9. [PSV] block coefficients against [N2](28) at finite g.
   10. Complex SYK: whether the symmetric channel is a single-scalar Robin problem.

Usage: python3 scripts/check_literature_inputs.py --dps 25
Exits non-zero if a check fails. Checks recorded as KNOWN-MISMATCH are printed differences
between sources, and are expected.
"""

import argparse
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_robin_identities import N2, kc_syk, response  # noqa: E402

RESULTS = []


def record(label, value, reference, tol=None):
    tol = tol if tol is not None else mp.mpf(10) ** (-mp.mp.dps + 8)
    diff = abs(mp.mpc(value) - mp.mpc(reference))
    scale = max(1, abs(mp.mpc(reference)))
    ok = diff / scale < tol
    RESULTS.append(ok)
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {label}: {mp.nstr(mp.mpc(value), 12)} vs {mp.nstr(mp.mpc(reference), 12)}")
    return ok


def note(label, text):
    print(f"  [KNOWN-MISMATCH] {label}: {text}")


# --------------------------------------------------------------------------------------
# Transcriptions of [PSV]
# --------------------------------------------------------------------------------------


class PSV:
    def __init__(self, q):
        self.q = mp.mpf(q)
        self.alpha = (self.q - 1) / self.q * mp.tan(mp.pi / (2 * self.q)) / (2 * mp.pi)  # (2.14)

    def ksd(self, h):  # (2.20)
        q, al = self.q, self.alpha
        return (
            al * mp.pi**2 * mp.cos(mp.pi * ((q + 2) / (4 * q) - h / 2)) * mp.gamma((q + 2) / (2 * q) - h)
        ) / (
            mp.sin(mp.pi / (2 * q)) * mp.gamma(1 / q) * mp.cos(mp.pi * (1 + q) / (2 * q)) * mp.gamma((1 + q) / q)
            * mp.sin(mp.pi * ((3 * q - 2) / (4 * q) - h / 2)) * mp.gamma((3 * q - 2) / (2 * q) - h)
        )

    def kad(self, h):  # (2.21)
        q, al = self.q, self.alpha
        return (
            al * mp.pi**2 * mp.sin(mp.pi * ((q + 2) / (4 * q) - h / 2)) * mp.gamma((q + 2) / (2 * q) - h)
        ) / (
            mp.sin(mp.pi / (2 * q)) * mp.gamma(1 / q) * mp.cos(mp.pi * (1 + q) / (2 * q)) * mp.gamma((1 + q) / q)
            * mp.cos(mp.pi * ((3 * q - 2) / (4 * q) - h / 2)) * mp.gamma((3 * q - 2) / (2 * q) - h)
        )

    def ka_matrix(self, h):  # (2.27)-(2.29)
        q, al = self.q, self.alpha
        common = (
            mp.pi**2 * mp.sin(mp.pi * (1 / (2 * q) - h / 2)) * mp.gamma(1 / q - h)
        ) / (
            mp.sin(mp.pi / (2 * q)) ** 2 * mp.gamma(1 / q) ** 2
            * mp.sin(mp.pi * ((2 * q - 1) / (2 * q) - h / 2)) * mp.gamma((2 * q - 1) / q - h)
        )
        k21 = (al / q) * (
            mp.pi**2 * mp.cos(mp.pi * ((q + 1) / (2 * q) - h / 2)) * mp.gamma((q + 1) / q - h)
        ) / (
            mp.cos(mp.pi * (q + 1) / (2 * q)) ** 2 * mp.gamma((q + 1) / q) ** 2
            * mp.cos(mp.pi * ((q - 1) / (2 * q) - h / 2)) * mp.gamma((q - 1) / q - h)
        )
        return mp.matrix([[-(q - 2) * al * common, -q * al * common], [k21, 0]])

    def ks_matrix(self, h):  # (2.33)-(2.35)
        q, al = self.q, self.alpha
        common = (
            mp.pi**2 * mp.cos(mp.pi * (1 / (2 * q) - h / 2)) * mp.gamma(1 / q - h)
        ) / (
            mp.sin(mp.pi / (2 * q)) ** 2 * mp.gamma(1 / q) ** 2
            * mp.cos(mp.pi * ((2 * q - 1) / (2 * q) - h / 2)) * mp.gamma((2 * q - 1) / q - h)
        )
        k21 = -(al / q) * (
            mp.pi**2 * mp.sin(mp.pi * ((q + 1) / (2 * q) - h / 2)) * mp.gamma((q + 1) / q - h)
        ) / (
            mp.cos(mp.pi * (q + 1) / (2 * q)) ** 2 * mp.gamma((q + 1) / q) ** 2
            * mp.sin(mp.pi * ((q - 1) / (2 * q) - h / 2)) * mp.gamma((q - 1) / q - h)
        )
        return mp.matrix([[(q - 2) * al * common, q * al * common], [k21, 0]])

    def ka_pm(self, h, sign):  # (2.30); sign=+1 is k^{a,+}
        q = self.q
        s = mp.mpf(sign)
        return (
            -s * mp.gamma(2 - 1 / q) * mp.gamma(1 - h / 2 - 1 / (2 * q)) * mp.gamma(1 / (2 * q) + h / 2)
            * mp.gamma(mp.mpf(1) / 2 - h + 1 / q - s / 2)
        ) / (
            mp.gamma(1 + 1 / q) * mp.gamma(1 + h / 2 - 1 / (2 * q)) * mp.gamma(1 / (2 * q) - h / 2)
            * mp.gamma(mp.mpf(3) / 2 - h - 1 / q - s / 2)
        )

    def ks_pm(self, h, sign):  # (2.36); sign=+1 is k^{s,+}
        q = self.q
        s = mp.mpf(sign)
        half = mp.mpf(1) / 2
        return (
            -s * mp.gamma(2 - 1 / q) * mp.gamma(half - h / 2 - 1 / (2 * q)) * mp.gamma(half + h / 2 + 1 / (2 * q))
            * mp.gamma(half - h + 1 / q - s / 2)
        ) / (
            mp.gamma(1 + 1 / q) * mp.gamma(half + h / 2 - 1 / (2 * q)) * mp.gamma(half - h / 2 + 1 / (2 * q))
            * mp.gamma(3 * half - h - 1 / q - s / 2)
        )

    def k0_psi_s(self, h):  # (4.10)
        q = self.q
        return -(mp.pi * mp.cos(mp.pi / (2 * q)) ** 2 * mp.gamma((q - 1) / q) ** 2) / (
            (mp.sin(mp.pi * h) + mp.sin(mp.pi / q)) * mp.gamma(-h - 1 / q + 2) * mp.gamma(h - 1 / q + 1)
        )

    def k0_psi_a(self, h):  # (4.11)
        q = self.q
        return (mp.pi * mp.cos(mp.pi / (2 * q)) ** 2 * mp.gamma((q - 1) / q) ** 2) / (
            (mp.sin(mp.pi / q) - mp.sin(mp.pi * h)) * mp.gamma(-h - 1 / q + 2) * mp.gamma(h - 1 / q + 1)
        )


def fgms_ks(qhat, h):  # [FGMS](6.2), first line
    qh = mp.mpf(qhat)
    x = 1 / (2 * qh)
    quarter = mp.mpf(1) / 4
    return (
        -mp.power(2, -1 + 2 / qh) * mp.gamma(2 - 1 / qh) / mp.gamma(1 + 1 / qh)
        * mp.gamma(quarter + x - h / 2) * mp.gamma(quarter + x + h / 2)
        / (mp.gamma(5 * quarter - x - h / 2) * mp.gamma(quarter - x + h / 2))
    )


def fgms_ka(qhat, h):  # [FGMS](6.2), second line
    qh = mp.mpf(qhat)
    x = 1 / (2 * qh)
    quarter = mp.mpf(1) / 4
    return (
        -mp.power(2, -1 + 2 / qh) * mp.gamma(2 - 1 / qh) / mp.gamma(1 + 1 / qh)
        * mp.gamma(3 * quarter + x - h / 2) * mp.gamma(-quarter + x + h / 2)
        / (mp.gamma(3 * quarter - x - h / 2) * mp.gamma(3 * quarter - x + h / 2))
    )


# --------------------------------------------------------------------------------------
# Transcriptions of [B18]
# --------------------------------------------------------------------------------------


def b18_kA(delta, h):  # (4.48), literal; 2*Delta = delta
    return (
        -mp.gamma(-delta) * mp.gamma(2 - delta) * mp.gamma(delta - h) * mp.gamma(delta + h)
        * mp.sin(mp.pi * delta) * (mp.sin(mp.pi * delta) - mp.sin(mp.pi * h)) / mp.pi**2
    )


def b18_kS(delta, h):  # (4.49), literal
    return (
        -mp.gamma(-delta) * mp.gamma(2 - delta) * mp.gamma(delta - h) * mp.gamma(delta + h)
        * mp.sin(mp.pi * delta) * (mp.sin(mp.pi * delta) + mp.sin(mp.pi * h)) / mp.pi**2
    )


def b18_kA_n0(delta, h):  # (B.8)
    return (
        mp.gamma(-delta) / mp.gamma(delta - 2) * mp.gamma(delta - h) * mp.gamma(delta + h - 1)
        * (mp.sin(mp.pi * h) - mp.sin(mp.pi * delta)) / mp.pi
    )


def b18_kS_n0(delta, h):  # (B.9)
    return (
        mp.gamma(1 - delta) / mp.gamma(delta - 1) * mp.gamma(delta - h) * mp.gamma(delta + h - 1)
        * (mp.sin(mp.pi * h) + mp.sin(mp.pi * delta)) / mp.pi
    )


def F_block(h, z):  # [B18](B.2)
    return mp.gamma(h) ** 2 / mp.gamma(2 * h) * z**h * mp.hyp2f1(h, h, 2 * h, z)


def psi_A(h, z):  # [B18](B.4), chi < 1
    c = mp.cos(mp.pi * h)
    return (2 / c) * (mp.cos(mp.pi * h / 2) ** 2 * F_block(h, z) - mp.sin(mp.pi * h / 2) ** 2 * F_block(1 - h, z))


def psi_S(h, z):  # [B18](B.5), chi < 1
    c = mp.cos(mp.pi * h)
    return (2 / c) * (-mp.sin(mp.pi * h / 2) ** 2 * F_block(h, z) + mp.cos(mp.pi * h / 2) ** 2 * F_block(1 - h, z))


def psi_integral(h, z, kind):
    """Direct evaluation of [B18](4.26) for 0 < chi=z < 1 and real 0 < h < 1."""

    def weight(y):
        base = z**h * abs(1 - y) ** (h - 1) / (abs(y) ** h * abs(z - y) ** h)
        if kind == "A":
            return base
        # ss' term of (4.26): sgn(chi) sgn(y) sgn(chi-y) sgn(1-y)
        return mp.sign(y) * mp.sign(z - y) * mp.sign(1 - y) * base

    pts = [-mp.inf, -1, 0, z / 2, z, (z + 1) / 2, 1, 2, mp.inf]
    return mp.quad(weight, pts)


# --------------------------------------------------------------------------------------
# Part 1
# --------------------------------------------------------------------------------------


def check_kernels():
    print("\n1. Component kernels: [PSV], [FGMS] against the draft")
    for q in (3, 5, 7):
        m = N2(q)
        psv = PSV(q)
        print(f" q={q}")
        for h in (mp.mpf("0.33"), mp.mpf("1.27"), mp.mpf("2.61")):
            # [N2] below (16): charged kernels. Superscripts as in [PSV].
            record(f"k^(a,d)({h}) = kA(h-1/2)", psv.kad(h), m.kA(h - mp.mpf(1) / 2))
            record(f"k^(s,d)({h}) = kS(h-1/2)", psv.ksd(h), m.kS(h - mp.mpf(1) / 2))
            # [PSV](2.31), (2.37) with [N2](16): the branch labels.
            record(f"k^(a,+)({h}) = kA(h)", psv.ka_pm(h, +1), m.kA(h))
            record(f"k^(a,-)({h}) = kS(h-1)", psv.ka_pm(h, -1), m.kS(h - 1))
            record(f"k^(s,+)({h}) = kS(h)", psv.ks_pm(h, +1), m.kS(h))
            record(f"k^(s,-)({h}) = kA(h-1)", psv.ks_pm(h, -1), m.kA(h - 1))
            # [N2](9), entry by entry. p=1 is [PSV] 'a' and p=0 is [PSV] 's', since [PSV]
            # labels the sector by the symmetry of the fermion component.
            omega = h - mp.mpf(1) / 2
            for p, mat in ((1, psv.ka_matrix(h)), (0, psv.ks_matrix(h))):
                kbar = m.kbar(p, omega)
                record(f"p={p} k11({h})", mat[0, 0], kbar[0, 0])
                record(f"p={p} k12({h})", mat[0, 1], kbar[0, 1] / (2 * m.delta))
                record(f"p={p} k21({h})", mat[1, 0], kbar[1, 0] * (2 * m.delta))
            # [FGMS](6.2): which label corresponds to which is checked here, not assumed.
            record(f"[FGMS] k_s({h}) = kS(h-1/2)", fgms_ks(q, h), m.kS(h - mp.mpf(1) / 2))
            record(f"[FGMS] k_a({h}) = kA(h-1/2)", fgms_ka(q, h), m.kA(h - mp.mpf(1) / 2))
        # [N2](15) checks, now from the literature side.
        record("k^(a,d)(3/2) = 1", psv.kad(mp.mpf(3) / 2), 1)
        record("k^(s,d)(3/2) = 1", psv.ksd(mp.mpf(3) / 2 + mp.mpf(10) ** (-mp.mp.dps + 6)), 1, tol=mp.mpf(10) ** -6)


# --------------------------------------------------------------------------------------
# Part 2
# --------------------------------------------------------------------------------------


def check_spectra():
    print("\n2. Printed q=3 spectra of [PSV] (2.38)-(2.43) against roots of the draft's kernels")
    m = N2(3)
    half = mp.mpf(1) / 2
    table = [
        ("(2.38) kA(h)", lambda h: m.kA(h), ["3.0659", "5.09488", "7.11311", "9.12623"]),
        ("(2.39) kA(h-1/2)", lambda h: m.kA(h - half), ["3.5659", "5.59488", "7.61311", "9.62623"]),
        ("(2.40) kA(h-1)", lambda h: m.kA(h - 1), ["4.0659", "6.09488", "8.11311", "10.1262"]),
        ("(2.41) kS(h)", lambda h: m.kS(h), ["2.82114", "4.74091", "6.69332", "8.66092"]),
        ("(2.42) kS(h-1/2)", lambda h: m.kS(h - half), ["3.32114", "5.24091", "7.19332", "9.16092"]),
        ("(2.43) kS(h-1)", lambda h: m.kS(h - 1), ["3.82114", "5.74091", "7.69332", "9.66092"]),
    ]
    for label, fn, roots in table:
        for printed in roots:
            root = mp.findroot(lambda h: fn(h) - 1, mp.mpf(printed))
            digits = len(printed.split(".")[1])
            record(f"{label} root {printed}", root, mp.mpf(printed), tol=mp.mpf(10) ** (-digits + 1))


# --------------------------------------------------------------------------------------
# Part 3
# --------------------------------------------------------------------------------------


def check_b18_kernel():
    print("\n3. [B18] kernel and the sign the draft says is misprinted")
    for q in (3, 5, 7):
        m = N2(q)
        d = m.delta
        print(f" q={q}")
        for h in (mp.mpf("0.37"), mp.mpf("1.41"), mp.mpf("2.73")):
            # The draft's claim: (4.48)-(4.49) carry the opposite overall sign.
            record(f"(4.48) literal = -kA({h})", b18_kA(d, h), -m.kA(h))
            record(f"(4.49) literal = -kS({h})", b18_kS(d, h), -m.kS(h))
            # [B18]'s own relation to N=0 (4.51) with (B.8) supports the draft's sign.
            record(f"(4.51)+(B.8) = +kA({h})", (d + h - 1) / (d - 2) * b18_kA_n0(d, h), m.kA(h))
            record(f"(4.52)+(B.9) = +kA({h})", (d - h - 1) / d * b18_kS_n0(d, -h), m.kA(h))
            record(f"(4.50) kA(h)=kS(-h)", m.kA(h), m.kS(-h))
        # [B18] states k=1 at h=1 for the charge multiplet. Only the corrected sign does.
        eps = mp.mpf(10) ** (-mp.mp.dps + 6)
        record("kA(1) = 1, draft sign", m.kA(1 + eps), 1, tol=mp.mpf(10) ** -6)
        record("(4.48) literal at h=1 is -1", b18_kA(d, 1 + eps), -1, tol=mp.mpf(10) ** -6)
    note("(4.48), (4.49)", "overall sign opposite to (4.51)+(B.8), [PSV], and k(1)=1. Draft's claim confirmed.")


# --------------------------------------------------------------------------------------
# Part 4
# --------------------------------------------------------------------------------------


def _xi(h, z):
    """xi_h = h (Psi^A_h - Psi^S_{h+1}), the definition [B18](4.38), from (B.4), (B.5)."""
    return h * (psi_A(h, z) - psi_S(h + 1, z))


def _phi(h, z):  # (4.11), equivalently (4.14)
    return F_block(h, z) - F_block(h + 1, z)


def check_b18_blocks():
    print("\n4. [B18] Eq. (4.39) and the seed normalization")
    print(" (B.4), (B.5) against the integral definition (4.26)")
    for h, z in ((mp.mpf("0.37"), mp.mpf("0.3")), (mp.mpf("0.62"), mp.mpf("0.55"))):
        record(f"Psi^A h={h} z={z}", psi_integral(h, z, "A"), psi_A(h, z), tol=mp.mpf(10) ** -8)
        record(f"Psi^S h={h} z={z}", psi_integral(h, z, "S"), psi_S(h, z), tol=mp.mpf(10) ** -8)
    print(" (4.39) with the draft's minus sign before phi_{-h}")
    for h, z in ((mp.mpf("0.37"), mp.mpf("0.3")), (mp.mpf("1.3"), mp.mpf("0.6"))):
        sec = 1 / mp.cos(mp.pi * h)
        corrected = h * ((1 + sec) * _phi(h, z) - (1 - sec) * _phi(-h, z))
        literal = h * ((1 + sec) * _phi(h, z) + (1 - sec) * _phi(-h, z))
        record(f"xi_h h={h} z={z}, minus sign", _xi(h, z), corrected)
        print(f"    literal (4.39) would give {mp.nstr(literal, 10)}; definition gives {mp.nstr(_xi(h, z), 10)}")
    note("(4.39)", "the printed plus sign disagrees with (4.38)+(B.4),(B.5). Draft's claim confirmed.")

    print(" Seed: leading coefficient of [B18](4.92) at the pole h=delta, with F0 -> z^delta")
    for q in (3, 5):
        m = N2(q)
        d = m.delta
        alpha = 4 * mp.pi * d / (m.nu * mp.tan(mp.pi * d / 2))  # (4.89) with 2 Delta = delta
        record(f"q={q} alpha/(4 pi^2) = D_delta", alpha / (4 * mp.pi**2), m.D)
        record(f"q={q} alpha/(8 pi) = N_delta", alpha / (8 * mp.pi), m.norm)

        def residue_at(c, h0):
            # upward-line weight from (4.92) with xi_h reduced by h -> -h as in the note;
            # numerator kernel uses the literal (4.48)-(4.49), as (4.88) is computed in App. G.
            def r(h):
                return (
                    -c * alpha / (4 * mp.pi)
                    * (mp.cot(mp.pi * h / 2) * b18_kA(d, h) - mp.tan(mp.pi * h / 2) * b18_kS(d, h))
                    * mp.gamma(h) ** 2 / mp.gamma(2 * h)
                )

            rad = d / 4
            f = lambda t: r(h0 + rad * mp.expjpi(2 * t)) * rad * mp.expjpi(2 * t) * 2j * mp.pi
            return mp.quad(f, [0, 0.5, 1]) / (2j * mp.pi)

        # Closing right gives -Res; the seed z^delta must have coefficient 1.
        record(f"q={q} factor 1/2 of (4.88) gives coefficient 1", -residue_at(mp.mpf(1) / 2, d), 1, tol=mp.mpf(10) ** -10)
        print(f"    without the 1/2, as printed in (4.90)-(4.92): {mp.nstr(-residue_at(mp.mpf(1), d), 10)}")
    note("(4.90)-(4.92), (4.106)", "prefactor alpha is 2x what (4.88) implies and what the seed needs.")


# --------------------------------------------------------------------------------------
# Part 5
# --------------------------------------------------------------------------------------


def ms_kc(q, h):  # [MS] v1 (3.73)
    q = mp.mpf(q)
    half = mp.mpf(1) / 2
    return (
        -(q - 1) * mp.gamma(3 * half - 1 / q) * mp.gamma(1 - 1 / q) / (mp.gamma(half + 1 / q) * mp.gamma(1 / q))
        * mp.gamma(1 / q + h / 2) / mp.gamma(3 * half - 1 / q - h / 2)
        * mp.gamma(half + 1 / q - h / 2) / mp.gamma(1 - 1 / q + h / 2)
    )


def check_ms():
    print("\n5. [MS] kernel against kc_syk, which was written from memory")
    for q in (4, 6, 10):
        for h in (mp.mpf("0.5") + 0.7j, mp.mpf("2.4"), mp.mpf("3.3")):
            record(f"q={q} h={mp.nstr(h, 3)}", kc_syk(q, h), ms_kc(q, h))
    h = mp.mpf("3.3")
    w = h - mp.mpf(1) / 2
    record("q=4 (3.75)", ms_kc(4, h), -mp.mpf(3) / 2 * mp.tan(mp.pi * w / 2) / w)
    record("k_c(2) = 1, q=4", ms_kc(4, mp.mpf(2)), 1)
    record("k_c(0) = -(q-1), q=6", ms_kc(6, mp.mpf(0)), -5)


# --------------------------------------------------------------------------------------
# Part 6
# --------------------------------------------------------------------------------------


def check_superfield_components():
    """A superprimary eigenfunction at tau_0 -> infinity is a function of <12> alone.

    Then f(<12>) = f(tau_12) - 2 theta-bar_1 theta_2 f'(tau_12) + ... ties the two
    components: the boson entry is proportional to the derivative of the fermion entry,
    so to (h - delta) for f ~ |tau|^(h-delta). The same must hold for the right eigenvector
    of the [PSV] 2x2 kernel, which knows nothing about superspace. See docs/derivations.md
    Sec. 5.

    The constant is fixed at h=0, where the eigenfunction is the propagator itself:
    (G_psi, G_b) ~ (sgn|t|^-delta, delta |t|^(-1-delta)), so c(0)=delta and c = -(h-delta).
    That is G_b = -d G_psi in the normalization of [PSV](2.10). Bulycheva's <12> has the
    coefficient 2, which corresponds to G_b twice as large; put into the Schwinger-Dyson
    equation this gives her J b^q = tan(pi delta/2)/(4 pi) instead of 1/(2 pi), Eq. (3.9).
    """
    print("\n6. Superfield structure of the top-branch eigenvector of the [PSV] matrix kernel")
    for q in (3, 5, 7):
        m = N2(q)
        psv = PSV(q)
        print(f" q={q}")
        for sector, mat_fn, branch in (("a", psv.ka_matrix, m.kA), ("s", psv.ks_matrix, m.kS)):
            for h in (mp.mpf("0.41"), mp.mpf("1.37"), mp.mpf("2.29")):
                mat = mat_fn(h)
                k = branch(h)
                # right eigenvector (1, c): k21 = k c
                c = mat[1, 0] / k
                record(f"{sector} h={h} eigenvector ratio c/(h-delta)", c / (h - m.delta), -1)
                # and it is an eigenvector: the first row closes too
                record(f"{sector} h={h} first row", mat[0, 0] + mat[0, 1] * c, k)


# --------------------------------------------------------------------------------------
# Part 7
# --------------------------------------------------------------------------------------


def _res(f, h0, rad):
    fn = lambda t: f(h0 + rad * mp.expjpi(2 * t)) * rad * mp.expjpi(2 * t)
    return mp.quad(fn, [0, 0.25, 0.5, 0.75, 1])


def _roots_G(m, p, g, omega_lo, omega_hi):
    def det(w):
        kb = m.kbar(p, w)
        return (1 - g * kb[0, 0]) * (1 - g * kb[1, 1]) - g**2 * kb[0, 1] * kb[1, 0]

    found = []
    grid = [omega_lo + (omega_hi - omega_lo) * i / 400 for i in range(401)]
    for a, b in zip(grid, grid[1:]):
        try:
            fa, fb = det(a), det(b)
        except ZeroDivisionError:
            continue
        if mp.im(fa) == 0 and mp.re(fa) * mp.re(fb) < 0 and abs(fa) < 50 and abs(fb) < 50:
            try:
                r = mp.findroot(det, (a, b), solver="anderson")
            except Exception:
                continue
            if abs(det(r)) < mp.mpf(10) ** -8 and not any(abs(r - s) < 1e-6 for s in found):
                found.append(mp.re(r))
    return found


def check_finite_g():
    """[N2](12) at finite g from the literature's own four-point function.

    [PSV](4.26): f^s = f0^s / [(1 - k^{a,+})(1 - k^{a,-})], with f0^s = 2 mu(h) k0^{psi,s} Phi^s
    where 2 mu(h) = (2h-1)/(pi tan(pi h/2)). Insert g by 1/(1 - g k) (the geometric series, as
    [T] does in (C.10)). Closing the contour gives block coefficients. [T](C.8) then gives
    P[2 N(h) B_h] = 2N(h) I(h) e^{-omega|tau|} = (1/mu(h)) e^{-omega|tau|}, so the factor
    tan(pi h/2) cancels:

        (P F_{p=1})(tau) = -2 sum_h Res[ k0^{psi,s}(h) / D^a(h) ] e^{-(h-1/2)|tau|}.

    The even sector is the same with k0^{psi,a}, D^s, I_a(h) from part 8, and overall sign +2.
    This is finite at h=1 even though the block coefficient vanishes there and I(1) diverges.
    [N2](12), with endpoint = mode/2, predicts 2 b^2 (P F)(tau) = sigma l_psi w_11(tau), with
    w = -Res G the coefficient of e^{-omega|tau|}.
    """
    print("\n7. Reduced four-point function at finite g: [PSV]+[T](C.8) against [N2](12)")
    half = mp.mpf(1) / 2
    for q in (3, 5):
        m = N2(q)
        psv = PSV(q)
        # J b^q = tan(pi delta/2)/(2 pi), [PSV](2.10); set J = 1, which [N2](12) is covariant under.
        b = mp.power(mp.tan(mp.pi * m.delta / 2) / (2 * mp.pi), m.delta)
        ell = 2 * mp.power(4, -m.nu) / b ** (q - 2)
        for g in (mp.mpf("0.5"), mp.mpf("0.9"), mp.mpf("0.99"), mp.mpf(1)):
            for p in (1, 0):
                if p == 0 and g < mp.mpf("0.9"):
                    continue  # below or near the even-sector threshold g_*; see derivations 3.3
                print(f" q={q} g={g}, p={p}")
                sigma = m.sigma(p)

                if p == 1:
                    D = lambda h: (1 - g * m.kA(h)) * (1 - g * m.kS(h - 1))
                    k0 = psv.k0_psi_s
                else:
                    D = lambda h: (1 - g * m.kS(h)) * (1 - g * m.kA(h - 1))
                    k0 = psv.k0_psi_a

                def G11(w):
                    return mp.inverse(mp.inverse(m.Dmat(p, w)) - g * m.Mmat(p))[0, 0]

                for w0 in _roots_G(m, p, g, mp.mpf("0.05"), mp.mpf(4)):
                    h0 = w0 + half
                    rad = mp.mpf("0.004")
                    # The sign is that of I(h) in each sector: see part 8.
                    res_psv = (-2 if p == 1 else 2) * _res(lambda h: k0(h) / D(h), h0, rad)
                    weight = -_res(G11, w0, rad)
                    record(
                        f"omega={mp.nstr(w0, 8)} (h={mp.nstr(h0, 8)})",
                        2 * b**2 * res_psv,
                        sigma * ell * weight,
                        tol=mp.mpf(10) ** -7,
                    )


def _block_a(h, chi):
    """The even-sector analogue of [T]'s B_h, from [PSV](3.7), (3.11): Phi^a split like [T](C.4)."""
    half = mp.mpf(1) / 2
    bh = lambda z: z**h * mp.hyp2f1(h, h, 2 * h, z)
    if chi < 0:
        return -bh(chi / (chi - 1))
    if chi < 1:
        return bh(chi)
    x = (chi - 2) / chi
    phi_a = (
        -2 * mp.gamma(1 - h / 2) * mp.gamma(h / 2 + half) / mp.sqrt(mp.pi)
        * x * mp.hyp2f1((2 - h) / 2, (h + 1) / 2, 3 * half, x**2)
    )
    return phi_a / (mp.gamma(h) ** 2 / mp.gamma(2 * h))


def _block_s(h, chi):
    """[T](C.4) with Psi_h from [PSV](3.6), as a control on the integration below."""
    half = mp.mpf(1) / 2
    bh = lambda z: z**h * mp.hyp2f1(h, h, 2 * h, z)
    if chi < 0:
        return bh(chi / (chi - 1))
    if chi < 1:
        return bh(chi)
    x = (chi - 2) / chi
    phi_s = (
        mp.gamma(half - h / 2) * mp.gamma(h / 2) / mp.sqrt(mp.pi)
        * mp.hyp2f1(h / 2, half - h / 2, half, x**2)
    )
    return phi_s / (mp.gamma(h) ** 2 / mp.gamma(2 * h))


def _centre_integral(block, h, tau):
    """(P B)(tau) of [T](C.8): integrate over u with chi = 1/(u^2 - sinh^2(tau/2)).

    Regions: chi < 0 for |u| < s, chi > 1 for s < |u| < c, 0 < chi < 1 for |u| > c.
    The singularities at |u| = s, c are logarithmic, handled by mp.quad's endpoint rule.
    """
    s = mp.sinh(tau / 2)
    c = mp.cosh(tau / 2)
    f = lambda u: block(h, 1 / (u**2 - s**2))
    return 2 * (mp.quad(f, [0, s]) + mp.quad(f, [s, c]) + mp.quad(f, [c, 2 * c, 10 * c, mp.inf]))


def check_centre_reduction():
    print("\n8. Centre integral of the blocks, odd sector ([T](C.8)) and even sector (new)")
    half = mp.mpf(1) / 2
    saved = mp.mp.dps
    mp.mp.dps = 15  # quadrature with logarithmic endpoints; 15 digits is enough here
    for h, tau in ((mp.mpf("1.3"), mp.mpf("0.7")), (mp.mpf("1.7"), mp.mpf("1.5")), (mp.mpf("2.4"), mp.mpf("0.9"))):
        w = h - half
        N = mp.gamma(h) ** 2 / (2 * mp.gamma(2 * h))
        pred_s = mp.pi * mp.tan(mp.pi * h / 2) / (2 * w * N) * mp.exp(-w * tau)
        # I_a(h) = -1/(2 mu_a(h) N(h)) with mu_a = (h-1/2)/(pi tan(pi (h-1)/2)), i.e. cot(pi h/2)
        # in place of the tan(pi h/2) of [T](C.8). It vanishes at h=1 and has a pole at h=2.
        pred_a = mp.pi * mp.cot(mp.pi * h / 2) / (2 * w * N) * mp.exp(-w * tau)
        record(f"odd, control: h={h} tau={tau}", _centre_integral(_block_s, h, tau), pred_s, tol=mp.mpf(10) ** -5)
        record(f"even: h={h} tau={tau}", _centre_integral(_block_a, h, tau), pred_a, tol=mp.mpf(10) ** -5)
    # The exceptional A-branch modes: at odd h, Phi^s has a pole living only on chi > 1 (and at
    # even h, Phi^a). Its residue is a Legendre polynomial in x = (chi-2)/chi, supported on chi > 1.
    print(" Step-function families: P[theta(chi>1) P_n((chi-2)/chi)] = 2/(2n+1) e^{-(n+1/2)|tau|}")
    for n in range(0, 5):
        for tau in (mp.mpf("0.4"), mp.mpf("1.3")):
            sh, ch = mp.sinh(tau / 2), mp.cosh(tau / 2)
            f = lambda u: mp.legendre(n, 1 - 2 * (u**2 - sh**2))
            record(f"n={n} tau={tau}", 2 * mp.quad(f, [sh, ch]), 2 / mp.mpf(2 * n + 1) * mp.exp(-(n + half) * tau), tol=mp.mpf(10) ** -10)
    print(" Residues of the partial waves at their poles, on chi > 1")
    chi = mp.mpf("1.6")
    x = (chi - 2) / chi
    rad = mp.mpf("0.01")
    phi_s = lambda h: (
        mp.gamma(half - h / 2) * mp.gamma(h / 2) / mp.sqrt(mp.pi) * mp.hyp2f1(h / 2, half - h / 2, half, x**2)
    )
    phi_a = lambda h: (
        -2 * mp.gamma(1 - h / 2) * mp.gamma(h / 2 + half) / mp.sqrt(mp.pi) * x * mp.hyp2f1((2 - h) / 2, (h + 1) / 2, 3 * half, x**2)
    )
    record("Res_{h=1} Phi^s = -2 P_0", _res(phi_s, mp.mpf(1), rad), -2, tol=mp.mpf(10) ** -10)
    record("Res_{h=3} Phi^s = -2 P_2", _res(phi_s, mp.mpf(3), rad), -2 * mp.legendre(2, x), tol=mp.mpf(10) ** -10)
    record("Res_{h=2} Phi^a = +2 P_1", _res(phi_a, mp.mpf(2), rad), 2 * mp.legendre(1, x), tol=mp.mpf(10) ** -10)
    record("Res_{h=4} Phi^a = +2 P_3", _res(phi_a, mp.mpf(4), rad), 2 * mp.legendre(3, x), tol=mp.mpf(10) ** -10)
    mp.mp.dps = saved


def check_block_coefficients():
    """[PSV] block coefficients against [N2](28), at finite g, for 0 < z < 1 with chi = z/(z-1).

    Since chi < 0 here, Phi^s -> +Phi^s(z) and Phi^a -> -Phi^a(z) by [PSV](3.11). With [N2]'s
    statement that the fermion covariance is -F and phi_H = F_H - F_{H+1}:

        coefficient of F_H(z)     at a root of kA(H) = 1/g:  -Res[2mu k0s/D^a] = -c_A(H),
        coefficient of F_{H+1}(z) at the same root:          +Res[2mu_a k0a/D^s] = +c_A(H),
    and with the roles of the two sectors swapped for kS.
    """
    print("\n9. Block coefficients: [PSV] against [N2](28)")
    half = mp.mpf(1) / 2
    for q in (3, 5):
        m = N2(q)
        psv = PSV(q)
        for g in (mp.mpf("0.7"), mp.mpf(1)):
            print(f" q={q} g={g}")
            two_mu = lambda h: (2 * h - 1) / (mp.pi * mp.tan(mp.pi * h / 2))
            two_mu_a = lambda h: (2 * h - 1) / (mp.pi * mp.tan(mp.pi * (h - 1) / 2))
            Da = lambda h: (1 - g * m.kA(h)) * (1 - g * m.kS(h - 1))
            Ds = lambda h: (1 - g * m.kS(h)) * (1 - g * m.kA(h - 1))
            odd = lambda h: two_mu(h) * psv.k0_psi_s(h) / Da(h)
            even = lambda h: two_mu_a(h) * psv.k0_psi_a(h) / Ds(h)
            rad = mp.mpf("0.002")
            for branch, k_fn, k_der in (
                ("A", m.kA, lambda H: mp.diff(m.kA, H)),
                ("S", m.kS, lambda H: mp.diff(m.kS, H)),
            ):
                grid = [mp.mpf("1.5") + mp.mpf(i) * mp.mpf("4.5") / 901 for i in range(902)]
                brackets = []
                for a, b in zip(grid, grid[1:]):
                    fa, fb = k_fn(a) - 1 / g, k_fn(b) - 1 / g
                    if fa * fb < 0 and abs(fa) < 20 and abs(fb) < 20:
                        brackets.append((a, b))
                for a, b in brackets:
                    H = mp.findroot(lambda x: k_fn(x) - 1 / g, (a, b), solver="anderson")
                    if branch == "A":
                        c = m.norm * mp.cot(mp.pi * H / 2) / (g**2 * k_der(H))
                        record(f"c_A H={mp.nstr(H, 8)}: F_H", _res(odd, H, rad), c, tol=mp.mpf(10) ** -7)
                        record(f"c_A H={mp.nstr(H, 8)}: F_H+1", _res(even, H + 1, rad), c, tol=mp.mpf(10) ** -7)
                    else:
                        c = -m.norm * mp.tan(mp.pi * H / 2) / (g**2 * k_der(H))
                        record(f"c_S H={mp.nstr(H, 8)}: F_H", -_res(even, H, rad), c, tol=mp.mpf(10) ** -7)
                        record(f"c_S H={mp.nstr(H, 8)}: F_H+1", -_res(odd, H + 1, rad), c, tol=mp.mpf(10) ** -7)


def check_complex_syk():
    """Complex SYK at zero charge: is the symmetric channel a single-scalar Robin problem?

    [T] covers the antisymmetric channel, k^A R_1 = (q-1) 4^-nu. If k^S R_0 is constant too, the
    symmetric channel is Thong's J=0 sector. It then has a + J = 2 Delta - 1/2 < 0 for q > 4,
    the same failure as the N=2 even sector, and a mode at h=1 (the U(1) charge).
    """
    print("\n10. Complex SYK: [B18](B.9), [GKST](B.3), [B17](4.29) and the J=0 response")
    half = mp.mpf(1) / 2
    for q in (4, 6, 8, 12):
        d2 = mp.mpf(2) / q  # 2 Delta, with Delta = 1/q
        nu = 1 - d2
        print(f" q={q}")
        for h in (mp.mpf("1.37"), mp.mpf("2.21"), mp.mpf("0.5") + 0.8j):
            gkst = (
                mp.gamma(d2 - h) * mp.gamma(d2 + h - 1) / (mp.gamma(d2 - 1) * mp.gamma(d2))
                * (1 + mp.sin(mp.pi * h) / mp.sin(mp.pi * d2))
            )
            record(f"[B18](B.9) = [GKST](B.3) h={mp.nstr(h, 3)}", b18_kS_n0(d2, h), gkst)
            record(f"[B18](B.8) = k_c h={mp.nstr(h, 3)}", b18_kA_n0(d2, h), ms_kc(q, h))
            record(f"k^S R_0 = 4^-nu, h={mp.nstr(h, 3)}", gkst * response(nu, 0, h - half), mp.power(4, -nu))
        eps = mp.mpf(10) ** (-mp.mp.dps + 8)
        record("k^S(1) = 1", b18_kS_n0(d2, 1 + eps), 1, tol=mp.mpf(10) ** -6)
        print(f"    a + J = {mp.nstr(d2 - half, 5)}")
    h = mp.mpf("2.21")
    w = h - half
    record("[B17](4.29), q=4", b18_kS_n0(mp.mpf(1) / 2, h), half / mp.tan(mp.pi * w / 2) / w)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dps", type=int, default=25, help="mpmath working precision")
    args = parser.parse_args()
    mp.mp.dps = args.dps
    print(f"mpmath {mp.__version__}, dps={args.dps}")

    check_kernels()
    check_spectra()
    check_b18_kernel()
    check_b18_blocks()
    check_ms()
    check_superfield_components()
    check_finite_g()
    check_centre_reduction()
    check_block_coefficients()
    check_complex_syk()

    failed = RESULTS.count(False)
    print(f"\n{len(RESULTS) - failed}/{len(RESULTS)} checks passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
