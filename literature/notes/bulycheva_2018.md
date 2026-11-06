# Bulycheva 2018

**Citation.** K. Bulycheva, *N = 2 SYK model in the superspace formalism*, JHEP 04 (2018) 036, arXiv:1801.09006v2. PDF: `papers/bulycheva_2018.pdf`. Tag: **[B18]**.

**Read:** Secs. 1–4, the start of Sec. 5, Sec. 7, Apps. A–F, and App. G through (G.7). The rest of Sec. 5 (retarded kernel) and most of Sec. 6 (the 2d model) were skimmed, as the project doesn't depend on them yet.

## Main question

Do the N=2 four-point calculation in manifest superspace, expanding in eigenfunctions of the $su(1,1|1)$ Casimir.

## Physical system

Chiral superfields $\Psi=\psi(\tau+\theta\bar\theta)+\theta b$, with a random holomorphic superpotential, odd $\hat q$ (3.1), $\hat q=q$ in project notation. Also a 2d analogue.

## Definitions

- $\langle12\rangle=\tau_{12}-2\bar\theta_1\theta_2-\theta_1\bar\theta_1-\theta_2\bar\theta_2$ (2.12), chiral in point 1 and antichiral in 2. It is the unique such invariant.
- $G(1|2)=\langle\bar\Psi(1)\Psi(2)\rangle=b\,\mathrm{sgn}\langle12\rangle/|\langle12\rangle|^{2\Delta}$, $2\Delta\hat q=1$, $4\pi Jb^{\hat q}=\tan\pi\Delta$ (3.5)–(3.9).
- Cross ratio $\chi=\langle12\rangle\langle34\rangle/(\langle14\rangle\langle32\rangle)$ (4.3). There's only one, and no nilpotent invariant, so $\mathcal W=\mathcal W(\chi)$ is purely bosonic (4.7). **This $\chi$ is the draft's $z$.** It is not PSV's $\chi$.

## Main results

- **Casimir** (4.9), (A.17): $\mathcal C=\chi^2(1-\chi)\partial_\chi^2+\chi(1-\chi)\partial_\chi$, eigenvalue $h^2$, symmetric under $h\to-h$ (4.10). Eigenfunctions $\varphi_h=\chi^hB(h,h)\,{}_2F_1(h,h;1+2h;\chi)=F_h-F_{h+1}$ (4.11), (4.14).
- **Eigenfunctions** from shadows: $\Xi^{AA},\Xi^{SS}$ are $T$-odd, $\Xi^{AS},\Xi^{SA}$ are $T$-even, and $\xi_h\equiv\Xi^{AS}_h=h(\Psi^A_h-\Psi^S_{h+1})$ (4.38).
- **Kernel** (4.48), (4.49). $k^A(h)=k^S(-h)$ (4.50), and $k^A(h)=k^A(-h)$ at integers (4.91). **Not even in $h$**, unlike the $N=0$ kernels.
- **Charge multiplet**: $h=1$, $\mathcal Q=R+\theta\bar Q+\bar\theta Q+\theta\bar\theta T$, dimensions $(1,\tfrac32,\tfrac32,2)$ (4.55).
- **Inner product** (4.57): $d\chi/[\chi(1-\chi)]$, not positive-definite. Not a Hilbert space, though positive on $\xi_h$: (4.72), (4.81). The constant and $\log\chi$ are orthogonal to the whole basis, (4.83), (4.84).
- **Contour**: the double pole at $h=0$ must be passed on the right, which the seed selects; on the left gives $\log\chi$ (Sec. 4.5).
- **Four-point function** (4.106), (4.107). Discrete sum and contour integral, the integer poles cancelling. Lyapunov exponent maximal.

## Printed errors found

Checked in `scripts/check_literature_inputs.py` and `scripts/check_supercasimir.py`. These are our findings and not an erratum. **The first two were claimed by the draft and are confirmed; the rest are new.**

| Eq. | Problem | Evidence |
|---|---|---|
| (4.48), (4.49) | Overall sign | Disagrees with her own (4.51)+(B.8), with [PSV], and with her statement $k(1)=1$ (literal gives $-1$) |
| (4.39) | Sign before $\varphi_{-h}$ should be minus | Disagrees with (4.38)+(B.4), (B.5); those we confirmed by direct integration of (4.26) |
| (4.90)–(4.92), (4.106) | Prefactor $\alpha$ should be $\alpha/2$ | Disagrees with (4.85)+(4.88), and the seed's leading coefficient comes out 2. (4.99) is correct |
| (A.9) | Should read $-L_{-1}L_1$ | As printed the operator isn't central. Reordered, it is, and gives (A.10), (A.16) |
| (A.3), (A.5), (A.7) | $Q/2$ should be $Q$ | As printed the algebra doesn't close with (A.4) |
| (A.15) | $Q=0$ holds for the pair only | $\langle12\rangle^{-2\Delta}$ is invariant only for $Q_1=-2\Delta$, $Q_2=+2\Delta$ |

None affects her final structure. With the corrected signs, (4.106) reorganizes exactly into [N2](26) (see `docs/derivations.md` Sec. 5.3).

**Normalization difference with [PSV].** Her $Jb^q=\tan(\pi\Delta)/(4\pi)$, against $1/(2\pi)$ in [PSV](2.10). This is a normalization of $b$: her $\langle12\rangle$ implies $G_b=-2\partial G_\psi$, whereas [PSV] has $G_b=-\partial G_\psi$. Putting a factor of 2 into the Schwinger–Dyson equation accounts for the difference (derived by hand). Kernels don't depend on this.

## Important equations for us

(2.12), (4.9), (4.11), (4.14), (4.38), (4.48)–(4.50), (4.55), (4.88), (4.106), (A.1)–(A.17).

## Numerical methods

Analytic; (4.99) is said to have been checked numerically.

## Limitations

The $h=1$ subspace is left to "the theory outside the conformal limit", with [Fu et al.] cited for it. Not resolved here.

## Relationship to our project

- **Superblocks and Casimir**: source of [N2](24), and of $\mathcal C_2$. Agree.
- **It answers B1.** Appendix A is enough to derive the N=2 relative cylinder (`docs/derivations.md` Sec. 6): at zero total supermomentum the pair Casimir is $L_0^2$, which is why the eigenvalue is $h^2$ and the shadow is $-h$. She does not make this statement; she computes the cross-ratio form directly.
- **The purely bosonic relative space**: counting $2|2$ pair coordinates minus $1|2$ supertranslations leaves $1|0$, so no odd relative coordinate. That is the same fact as her (4.7), one level down.
