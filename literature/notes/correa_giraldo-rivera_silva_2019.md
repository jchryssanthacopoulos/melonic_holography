# Correa, Giraldo-Rivera, Silva 2019

**Citation.** D. H. Correa, V. I. Giraldo-Rivera, G. A. Silva, *Supersymmetric mixed boundary conditions in AdS$_2$ and DCFT$_1$ marginal deformations*, JHEP 03 (2020) 010, arXiv:1910.04225v3. PDF: `papers/correa_giraldo-rivera_silva_2019.pdf`. Tag: **[CGS]**.

**Read:** Secs. 1, 2, 3.1, 3.2, 5, and App. A. Secs. 3.3 and 4 and App. B were not read.

## Main question

Which mixed boundary conditions for scalars and fermions on AdS$_2$ preserve supersymmetry, and which is dual to the ABJM interpolating Wilson loops?

## Physical system

Free fields on rigid AdS$_2$, as fluctuations of an open string: an N=2 multiplet of a complex scalar and a Dirac fermion (2.10).

## Important definitions

- **Multiplet condition** (2.12): $m_B^2=m_F^2-m_F$, with $\delta\psi=(-i\gamma^\alpha\partial_\alpha\phi+m_F\phi)\epsilon$, $\delta\phi=\bar\epsilon\psi$. Four real supercharges.
- Quantization windows: alternate quantization for $-\tfrac14\le m_B^2\le\tfrac34$ and $|m_F|\le\tfrac12$.
- **All boundary-condition analysis in the paper is for massless fields**, for which $\phi\sim\alpha+y\beta$ and $\psi\sim y^{1/2}(\alpha^\psi+\beta^\psi)$ with $P_\mp$ projections.

## Main analytical results

| Type | Boson | Fermion | Supercharges | Scale invariant |
|---|---|---|---|---|
| Dirichlet (2.20) | $\alpha=0$ | $\alpha^\psi=0$ | 4 | yes |
| Neumann (2.21) | $\beta=0$ | $\beta^\psi=0$ | 4 | yes |
| I (2.24) | $\chi\alpha-\beta=0$ | $\chi\alpha^\psi-\gamma_5\dot\beta^\psi=0$ | **2** ($\xi_0$ only) | **no** |
| II (2.31) | $i\chi\dot\alpha-\beta=0$ | $i\chi\alpha^\psi-\gamma_5\beta^\psi=0$ | 4 | yes |
| III (2.36)–(2.40) | II plus constraints at $t_0,\tilde t_0$ | same | 4 | yes |

Boundary terms: (2.41) for I, (2.48) for II, (2.52) for III. Fermions need none in these cases.

Pattern: **in the supersymmetric pair, one member has a derivative and the other doesn't.** A plain Robin condition on the scalar pairs with a derivative condition on the fermion, and a plain mixed condition on the fermion pairs with a derivative condition on the scalar.

Global AdS$_2$ (App. A.2): Killing spinor (A.7) $\epsilon\propto(\cos\tfrac t2-i\gamma_0\sin\tfrac t2)\xi_0$, so **the supercharges carry global energy $\pm\tfrac12$**. Type II works in global coordinates and with both boundaries (A.21)–(A.23).

## Numerical methods

None.

## Limitations

Massless fields only. Poincaré time for most of it. No SYK; their closing line suggests supersymmetric SYK as a possible application ("...as holographic double trace deformations in other setups... for example supersymmetric SYK models").

## Relationship to our project

The draft cites this as background on boundary data. It bears on C4 in three ways:

1. **Mass pattern (checked by hand).** $m_b^2=\mu(\mu-1)$ is (2.12) at $m_F=\mu$, and $m_\psi^2=\mu(\mu+1)$ is (2.12) at $m_F=-\mu$. Both draft scalars are in the window since $0<\mu<\tfrac12$. So the masses are compatible with AdS$_2$ N=2 multiplets, but **the draft's scalars are not obviously the superpartners in this sense**: the draft gets them by squaring a Dirac equation, (19), with no scalar-fermion action.
2. **Energy shift.** The Killing spinor's energy $\pm\tfrac12$ matches the component shift $\omega_{\rm comp}=H\mp\tfrac12$ relative to $\Omega=H$ (see `docs/derivations.md` Sec. 6). The draft's test [N2](35) uses the isometry $E_+$, which shifts by 1. **The analogous test with a Killing spinor has not been done.**
3. **Candidate structure.** The draft's Dirac condition $b_F=gσL_qa_F$ is algebraic, like the fermion in Type II. Through [N2](20) it becomes frequency-dependent on the scalars, reminiscent of $i\chi\dot\alpha-\beta=0$. This is only an analogy so far, since [CGS] has $m=0$, and is worth testing with $m_F=\mu$ (task in `docs/todo.md`).

Also note that Type I shows a plain Robin condition *can* keep half the supercharges. So [N2](35) does not settle how many survive.
