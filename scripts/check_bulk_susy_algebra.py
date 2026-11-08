#!/usr/bin/env python3
"""Supersymmetry of the free AdS2 multiplet (Phi_b, Psi, Phi_psi), symbolically.

Strip metric ds^2 = sec^2(theta) (-dt^2 + dtheta^2), mostly plus, with frames gamma^0 = i sigma_1,
gamma^1 = -sigma_3, so that (nabla-slash - mu) Psi = 0 with Psi = cos^(1/2)(theta) psi is exactly
i d_t psi = (J2 d_theta + mu sec(theta) sigma_1) psi, the draft's [N2] Dirac Hamiltonian.

Tests:
  1. epsilon_A, epsilon_B solve the Killing equation nabla_a eps = (1/2) gamma_a eps.
  2. (nabla-slash - mu)(d-slash + mu) Phi eps = (box - mu(mu-1)) Phi eps, and
     (nabla-slash - mu) gamma3 (d-slash - mu) Phi eps = -gamma3 (box - mu(mu+1)) Phi eps,
     for arbitrary Phi(t, theta). These are what make the free action invariant, for any beta in
       delta Psi = [(d-slash + mu) Phi_b + beta gamma3 (d-slash - mu) Phi_psi] eps.
  3. Closure on Phi_b: the Phi_psi cross term vanishes if and only if beta is imaginary, and the
     normalization on Phi_psi requires |beta| = 1. So beta = +i or -i, two su(1,1|1) embeddings.
  4. The result is an AdS2 Killing vector: d_t for (A, A) and (B, B), E_pm for (A, B) and (B, A).

Usage: python3 scripts/check_bulk_susy_algebra.py
"""

import sys

import sympy as sp

t, th = sp.symbols("t theta", real=True)
mu = sp.symbols("mu", positive=True)
I = sp.I
s1 = sp.Matrix([[0, 1], [1, 0]])
s2 = sp.Matrix([[0, -I], [I, 0]])
s3 = sp.Matrix([[1, 0], [0, -1]])
one = sp.eye(2)
g0 = I * s1  # gamma^0
g1 = -s3  # gamma^1
g3 = g0 * g1  # chirality, = -sigma_2 up to the factor checked below
Lam = -sp.log(sp.cos(th))  # metric = exp(2 Lam) eta

RESULTS = []


def check(label, ok):
    RESULTS.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")


def is_zero(expr):
    return sp.simplify(expr) == 0


def mat_zero(M):
    return all(is_zero(sp.simplify(x.rewrite(sp.cos))) for x in M)


def spin_cov(eps, a):
    """nabla_a eps in the conformal frame: d_a + (1/2) gamma_a^b d_b Lam, with a = 0, 1 curved."""
    x = (t, th)[a]
    dL = [sp.diff(Lam, t), sp.diff(Lam, th)]
    glow = [-g0, g1]  # gamma_a, flat index lowered with eta = diag(-1, 1)
    gup = [g0, g1]
    res = sp.diff(eps, x)
    for b in range(2):
        if b == a:
            continue
        gam = (glow[a] * gup[b] - gup[b] * glow[a]) / 2  # gamma_a^b
        res += sp.Rational(1, 2) * gam * eps * dL[b]
    return res


def killing_residual(eps, kappa):
    glow = [-g0, g1]
    return [sp.simplify(spin_cov(eps, a) - kappa * sp.exp(Lam) * glow[a] * eps) for a in range(2)]


def dslash(phi):
    """d-slash of a scalar, as a matrix: cos(theta) gamma^a d_a phi."""
    return sp.cos(th) * (g0 * sp.diff(phi, t) + g1 * sp.diff(phi, th))


def box(phi):
    return sp.cos(th) ** 2 * (-sp.diff(phi, t, 2) + sp.diff(phi, th, 2))


def dirac(Psi):
    """nabla-slash on a spinor: cos(theta) gamma^a (d_a + (1/2) d_a Lam)."""
    return sp.cos(th) * (g0 * sp.diff(Psi, t) + g1 * (sp.diff(Psi, th) + sp.Rational(1, 2) * sp.tan(th) * Psi))


# Killing spinors. w3, w2 are the real rotated vectors in the notes, times cos^(-1/2).
c_p = sp.cos(sp.pi / 4 + th / 2)
c_m = sp.cos(sp.pi / 4 - th / 2)
pref = sp.cos(th) ** sp.Rational(-1, 2)
eps_A = pref * sp.Matrix([c_p, c_m]) * sp.exp(-I * t / 2)
eps_B = pref * sp.Matrix([-c_p, c_m]) * sp.exp(I * t / 2)


def bar(eps):
    """Dirac conjugate: eps^dagger sigma_1. For real t, theta in the strip."""
    conj = eps.applyfunc(lambda x: x.subs(I, -I))
    return conj.T * s1


def main():
    print("gamma3 = gamma^0 gamma^1 =", g3.tolist())
    print("1. Killing equation")
    for name, eps in (("A", eps_A), ("B", eps_B)):
        good = None
        for kappa in (sp.Rational(1, 2), -sp.Rational(1, 2)):
            res = killing_residual(eps, kappa)
            if all(all(sp.simplify(sp.trigsimp(x)) == 0 for x in r) for r in res):
                good = kappa
        print(f"    eps_{name} solves nabla eps = kappa gamma eps with kappa = {good}")
        check(f"eps_{name} is Killing with kappa = 1/2", good == sp.Rational(1, 2))

    print("2. Intertwining relations for arbitrary Phi(t, theta)")
    Phi = sp.Function("Phi")(t, th)
    for name, eps in (("A", eps_A), ("B", eps_B)):
        chi1 = (dslash(Phi) + mu * Phi * one) * eps
        lhs1 = dirac(chi1) - mu * chi1
        rhs1 = (box(Phi) - mu * (mu - 1) * Phi) * eps
        check(f"(nabla - mu)(d + mu) Phi eps_{name} = (box - mu(mu-1)) Phi eps_{name}",
              all(sp.simplify(sp.trigsimp(x.doit())) == 0 for x in (lhs1 - rhs1)))
        chi2 = g3 * (dslash(Phi) - mu * Phi * one) * eps
        lhs2 = dirac(chi2) - mu * chi2
        rhs2 = -g3 * (box(Phi) - mu * (mu + 1) * Phi) * eps
        check(f"(nabla - mu) g3 (d - mu) Phi eps_{name} = -g3 (box - mu(mu+1)) Phi eps_{name}",
              all(sp.simplify(sp.trigsimp(x.doit())) == 0 for x in (lhs2 - rhs2)))

    print("3. Closure: delta Psi = -M eps, delta Phi_b = Psibar eps + ebar Psi,")
    print("   delta Phi_psi = -beta Psibar g3 eps + conj(beta) ebar g3 Psi, as the action requires.")
    print("   [delta_1, delta_2] Phi_b = ebar_1 (M - Mbar) eps_2 - (1 <-> 2), Mbar = s1 M^dagger s1.")
    Pb = sp.Function("Phi_b")(t, th)
    Pp = sp.Function("Phi_psi")(t, th)

    def Mop(b):
        return dslash(Pb) + mu * Pb * one + b * g3 * (dslash(Pp) - mu * Pp * one)

    def Mbar(b_conj):
        # conjugating by hand, with real fields: i -> -i in dslash, beta -> conj(beta)
        dag = lambda phi: sp.cos(th) * (-I * s1 * sp.diff(phi, t) - s3 * (-1) * (-1) * sp.diff(phi, th))
        Mdag = dag(Pb) + mu * Pb * one + b_conj * (dag(Pp) - mu * Pp * one) * g3
        return s1 * Mdag * s1

    for b_val, b_conj, tag, closes in ((I, -I, "beta=+i", True), (-I, I, "beta=-i", True), (sp.Integer(1), sp.Integer(1), "beta=1", False)):
        diff = (Mop(b_val) - Mbar(b_conj)).applyfunc(sp.simplify)
        pure = (2 * dslash(Pb) - diff).applyfunc(sp.simplify)
        check(f"{tag}: M - Mbar = 2 d-slash Phi_b, no Phi_psi {'(closes)' if closes else 'is FALSE (does not close)'}",
              all(x == 0 for x in pure) == closes)

    print("4. The Killing vectors xi^a = ebar_1 cos(theta) gamma^a eps_2")
    def xi(e1, e2):
        return [sp.simplify(sp.trigsimp((bar(e1) * (sp.cos(th) * g) * e2)[0])) for g in (g0, g1)]

    table = {"(A,A)": xi(eps_A, eps_A), "(B,B)": xi(eps_B, eps_B), "(A,B)": xi(eps_A, eps_B), "(B,A)": xi(eps_B, eps_A)}
    for k, v in table.items():
        print(f"    {k}: xi^t = {v[0]}, xi^theta = {v[1]}")
    check("(A,A) and (B,B) generate d_t", table["(A,A)"][1] == 0 and table["(B,B)"][1] == 0 and table["(A,A)"][0].has(t) is False)
    E = table["(A,B)"]
    conf = [sp.simplify(sp.diff(E[0], t) - sp.diff(E[1], th)), sp.simplify(sp.diff(E[0], th) - sp.diff(E[1], t))]
    check("(A,B) is a conformal Killing vector of the strip, so an isometry of AdS2 with E^theta = 0 at the boundary",
          all(c == 0 for c in conf) and sp.simplify(E[1].subs(th, sp.pi / 2)) == 0)
    check("(A,B) raises energy by one unit: xi ~ exp(i t)", sp.simplify(sp.diff(E[0], t) / E[0]) in (I, -I))

    failed = RESULTS.count(False)
    print(f"\n{len(RESULTS) - failed}/{len(RESULTS)} checks passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
