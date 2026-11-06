#!/usr/bin/env python3
"""Reduce the su(1,1|1) pair Casimir at zero total supermomentum.

Implements the generators of Bulycheva, arXiv:1801.09006v2, Eqs. (A.1)-(A.8) as differential
operators on a Grassmann algebra, with coefficients handled by SymPy. Residuals are evaluated
at a generic rational point (q=3 values), not simplified symbolically. Tests:

    1. Printed Eq. (A.9) has -L_1 L_{-1}, translation on the left. Which ordering is central?
    2. The R-charges: (A.15) sets Q=0, but only the total vanishes. Closure of the algebra
       and invariance of the two-point function fix Q/2 -> Q in (A.3), (A.5), (A.7), and
       Q_1 = -2 Delta, Q_2 = +2 Delta.
    3. Functions of <12> are annihilated by the total L_1, G_{-1/2}, Gbar_{-1/2}, J_0.
    4. On them the pair Casimir equals L_0^2 = (x d_x + 2 Delta)^2, with x = <12>.
    5. The component content of f(<12>).
    6. The cross-ratio form (A.16)-(A.17), testing the implementation against the paper.

Usage: python3 scripts/check_supercasimir.py
Exits non-zero if a check fails. See docs/derivations.md Sec. 5.
"""

import sys
from itertools import product

import sympy as sp

ODD = ("t1", "tb1", "t2", "tb2", "t3", "tb3", "t4", "tb4")
INDEX = {name: i for i, name in enumerate(ODD)}


class G:
    """Element of the Grassmann algebra: {sorted tuple of odd names: SymPy coefficient}."""

    def __init__(self, terms=None):
        self.terms = {k: v for k, v in (terms or {}).items() if v != 0}

    def __add__(self, other):
        other = lift(other)
        out = dict(self.terms)
        for k, v in other.terms.items():
            out[k] = out.get(k, 0) + v
        return G(out)

    __radd__ = __add__

    def __neg__(self):
        return G({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + (-lift(other))

    def __mul__(self, other):
        other = lift(other)
        out = {}
        for (a, x), (b, y) in product(self.terms.items(), other.terms.items()):
            if set(a) & set(b):
                continue
            seq = list(a) + list(b)
            sign = 1
            for i in range(len(seq)):
                for j in range(i + 1, len(seq)):
                    if INDEX[seq[i]] > INDEX[seq[j]]:
                        sign = -sign
            key = tuple(sorted(seq, key=INDEX.get))
            out[key] = out.get(key, 0) + sign * x * y
        return G(out)

    def __rmul__(self, other):
        return lift(other) * self

    def expand(self):
        return G({k: sp.simplify(v) for k, v in self.terms.items()})

    def is_zero(self):
        return all(sp.simplify(v) == 0 for v in self.terms.values())


def lift(x):
    return x if isinstance(x, G) else G({(): sp.sympify(x)})


def odd(name):
    return G({(name,): sp.Integer(1)})


def d_odd(f, name):
    """Left derivative with respect to an odd variable."""
    out = {}
    for key, val in f.terms.items():
        if name in key:
            i = key.index(name)
            rest = key[:i] + key[i + 1 :]
            out[rest] = out.get(rest, 0) + (-1) ** i * val
    return G(out)


def d_even(f, var):
    return G({k: sp.diff(v, var) for k, v in f.terms.items()})


def generators(n, delta, Q=0, qfac=1):
    """Eqs. (A.1)-(A.8), for the particle with even coordinate tau_n and R-charge Q.

    As printed, (A.3), (A.5) and (A.7) have Q/2, which is qfac=1/2. That does not close with
    (A.4); qfac=1 does. See test 2 in main().
    """
    tau = sp.Symbol(f"tau{n}", real=True)
    th, tb = odd(f"t{n}"), odd(f"tb{n}")
    dth = lambda f: d_odd(f, f"t{n}")
    dtb = lambda f: d_odd(f, f"tb{n}")
    dt = lambda f: d_even(f, tau)
    half = sp.Rational(1, 2)
    Q = sp.sympify(Q)
    gens = {}
    gens["L0"] = lambda f: -tau * dt(f) - half * (th * dth(f)) - half * (tb * dtb(f)) - delta * f
    gens["L1"] = lambda f: -dt(f)
    gens["Lm1"] = lambda f: (
        -tau**2 * dt(f) - tau * (th * dth(f)) - tau * (tb * dtb(f)) - 2 * tau * delta * f - (Q * qfac) * (th * tb * f)
    )
    gens["J0"] = lambda f: -(th * dth(f)) + tb * dtb(f) + Q * f
    gens["Gp"] = lambda f: tau * dtb(f) - tau * (th * dt(f)) - (2 * delta + Q * qfac) * (th * f) - (th * tb * dtb(f))
    gens["Gm"] = lambda f: dtb(f) - th * dt(f)
    gens["Gbp"] = lambda f: tau * dth(f) - tau * (tb * dt(f)) - (2 * delta - Q * qfac) * (tb * f) + (th * tb * dth(f))
    gens["Gbm"] = lambda f: dth(f) - tb * dt(f)
    return gens


def casimir(gens, ordering):
    """(A.9). ordering='printed' is -L1 L-1 as printed; 'right' is -L-1 L1, translation acting first."""

    def apply(f):
        out = gens["L0"](gens["L0"](f))
        out = out - sp.Rational(1, 4) * gens["J0"](gens["J0"](f))
        if ordering == "printed":
            out = out - gens["L1"](gens["Lm1"](f))
        else:
            out = out - gens["Lm1"](gens["L1"](f))
        half = sp.Rational(1, 2)
        out = out + half * gens["Gp"](gens["Gbm"](f)) + half * gens["Gbp"](gens["Gm"](f))
        return out

    return apply


def total(name, particles):
    def apply(f):
        out = lift(0)
        for gens in particles:
            out = out + gens[name](f)
        return out

    return apply


def pair_casimir(particles):
    names = ("L0", "L1", "Lm1", "J0", "Gp", "Gm", "Gbp", "Gbm")
    return casimir({n: total(n, particles) for n in names}, "right")


RESULTS = []


def check(label, ok):
    RESULTS.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")


def compose(func, x):
    """Taylor expand func(x) about the body of x; func is a SymPy function of a Symbol."""
    body = x.terms.get((), 0)
    soul = G({k: v for k, v in x.terms.items() if k != ()})
    s = sp.Symbol("s_", positive=True)
    out = lift(0)
    power = lift(1)
    fact = 1
    for n in range(0, 9):
        coeff = sp.diff(func(s), s, n).subs(s, body) / fact
        out = out + lift(coeff) * power
        power = power * soul
        fact *= n + 1
    return out


POINT = {
    sp.Symbol("tau1", real=True): sp.Rational(23, 10),
    sp.Symbol("tau2", real=True): sp.Rational(1, 7),
    sp.Symbol("tau3", real=True): sp.Rational(-5, 3),
    sp.Symbol("tau4", real=True): sp.Rational(-41, 7),
}


def size(f):
    """Largest coefficient at a generic point. Numerical, since symbolic zero tests can stall."""
    return max((abs(sp.N(v.subs(POINT), 30)) for v in f.terms.values()), default=sp.Integer(0))


def angle(i, j):
    """<ij> of (2.12): tau_ij - 2 thetabar_i theta_j - theta_i thetabar_i - theta_j thetabar_j."""
    ti, tj = sp.Symbol(f"tau{i}", real=True), sp.Symbol(f"tau{j}", real=True)
    return lift(ti - tj) - 2 * (odd(f"tb{i}") * odd(f"t{j}")) - (odd(f"t{i}") * odd(f"tb{i}")) - (odd(f"t{j}") * odd(f"tb{j}"))


def test_superfunction():
    tau = sp.Symbol("tau1", real=True)
    t, tb = odd("t1"), odd("tb1")
    return (
        lift(tau ** sp.Rational(3, 5))
        + t * lift(tau**2 + 1)
        + tb * lift(sp.exp(-tau))
        + t * tb * lift(sp.log(tau))
    )


def main():
    D = sp.Rational(1, 6)  # q=3: 2 Delta q = 1
    tol = 1e-20
    names = ("L0", "L1", "Lm1", "J0", "Gp", "Gm", "Gbp", "Gbm")

    print("1. One-particle Casimir at Q=0. Which ordering of (A.9) is central?")
    gens = generators(1, D)
    f = test_superfunction()
    for ordering in ("printed", "right"):
        cas = casimir(gens, ordering)
        worst = max(size(cas(gens[n](f)) - gens[n](cas(f))) for n in names)
        print(f"    ordering '{ordering}': largest [C, generator] on a test superfunction = {sp.N(worst, 6)}")
        if ordering == "printed":
            check("printed order -L_1 L_{-1} is not central with these generators", worst > 1e-3)
        else:
            check("order -L_{-1} L_1, with translation acting first, is central", worst < tol)
    plain = lift(sp.Symbol("tau1", real=True) ** sp.Rational(3, 5))
    check("(A.10): C f(tau) = Delta^2 f(tau) at Q=0", size(casimir(gens, "right")(plain) - D**2 * plain) < tol)

    print("2. R-charges. (A.15) sets Q=0, but only the total can vanish: |Q_i| = 2 Delta")
    Q = 2 * D

    def closure(g):
        return max(
            size(g["Gp"](g["Gbm"](f)) + g["Gbm"](g["Gp"](f)) - 2 * g["L0"](f) + g["J0"](f)),
            size(g["Gbp"](g["Gm"](f)) + g["Gm"](g["Gbp"](f)) - 2 * g["L0"](f) - g["J0"](f)),
            size(g["Gp"](g["Gbp"](f)) + g["Gbp"](g["Gp"](f)) - 2 * g["Lm1"](f)),
            size(g["Gm"](g["Gbm"](f)) + g["Gbm"](g["Gm"](f)) - 2 * g["L1"](f)),
        )

    check("printed Q/2 in (A.3), (A.5), (A.7) does not close", closure(generators(1, D, Q, sp.Rational(1, 2))) > 1e-3)
    g_fixed = generators(1, D, Q, 1)
    check("Q/2 -> Q closes", closure(g_fixed) < tol)
    cas_q = casimir(g_fixed, "right")
    check("and the Casimir is central at Q != 0", max(size(cas_q(g_fixed[n](f)) - g_fixed[n](cas_q(f))) for n in names) < tol)
    x12 = angle(1, 2)
    two_pt = compose(lambda u: u ** (-2 * D), x12)
    pair = (generators(1, D, -Q), generators(2, D, Q))
    worst = max(size(total(n, pair)(two_pt)) for n in names)
    check("<12>^(-2 Delta) is invariant under all generators for Q_1 = -2 Delta, Q_2 = +2 Delta", worst < tol)
    wrong = (generators(1, D, Q), generators(2, D, -Q))
    check("but not for the opposite assignment", max(size(total(n, wrong)(two_pt)) for n in names) > 1e-3)

    print("3. Zero total supermomentum")
    h = sp.Rational(7, 5)
    fx = compose(lambda u: u ** (h - 2 * D), x12)
    for name in ("L1", "Gm", "Gbm", "J0"):
        check(f"total {name} annihilates f(<12>)", size(total(name, pair)(fx)) < tol)

    print("4. Reduced Casimir on f(<12>) = <12>^(h - 2 Delta)")
    cas = casimir({n: total(n, pair) for n in names}, "right")
    check("C_pair f(<12>) = h^2 f(<12>)", size(cas(fx) - h**2 * fx) < tol)
    l0 = total("L0", pair)
    check("C_pair f(<12>) = L_0^2 f(<12>)", size(cas(fx) - l0(l0(fx))) < tol)
    for test, label in ((lambda u: u ** sp.Rational(3, 7), "u^(3/7)"), (lambda u: sp.exp(-u), "exp(-u)"), (lambda u: sp.log(u), "log u")):
        lhs = l0(compose(test, x12))
        target = compose(lambda u: -(u * sp.diff(test(u), u) + 2 * D * test(u)), x12)
        check(f"L_0 f(<12>) = -(x f' + 2 Delta f)(<12>) for f = {label}", size(lhs - target) < tol)

    print("5. Component content: f(<12>) gives boson entry -2 f'(tau_12)")
    tau1, tau2 = sp.Symbol("tau1", real=True), sp.Symbol("tau2", real=True)
    expected = 2 * sp.diff((tau1 - tau2) ** (h - 2 * D), tau1)
    # key order follows ODD, which puts tb1 before t2
    tb1t2 = fx.terms.get(("tb1", "t2"), fx.terms.get(("t2", "tb1"), 0))
    print(f"    coefficient of thetabar_1 theta_2 = {sp.N(tb1t2.subs(POINT), 10)}, -2 f' = {sp.N((-expected).subs(POINT), 10)}")
    check("coefficient of thetabar_1 theta_2 in f(<12>) is -2 f'(tau_12)", abs(sp.N((tb1t2 + expected).subs(POINT), 30)) < tol)

    print("6. Cross-ratio form (A.16)-(A.17), as a test of the implementation against the paper")
    p = sp.Rational(3, 5)
    pref = compose(lambda u: u ** (-2 * D), x12)
    chi = angle(1, 2) * angle(3, 4) * compose(lambda u: 1 / u, angle(1, 4) * angle(3, 2))
    fchi = compose(lambda u: u**p, chi)
    full = pref * fchi
    target_f = lambda u: u**2 * (1 - u) * sp.diff(u**p, u, 2) + u * (1 - u) * sp.diff(u**p, u)
    target = pref * compose(target_f, chi)
    residual = cas(full) - target
    worst = size(residual)
    print(f"    largest coefficient of the residual at a generic point: {sp.N(worst, 6)}")
    check("(A.16)-(A.17) reproduced, including all odd terms", worst < tol)

    failed = RESULTS.count(False)
    print(f"\n{len(RESULTS) - failed}/{len(RESULTS)} checks passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
