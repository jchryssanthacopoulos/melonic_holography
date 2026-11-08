# Thong 2026: Melonic ladder holography with a single free bulk field

## Full citation

K. H. Thong, *Melonic ladder holography with a single free bulk field in AdS$_{d+1}$*, arXiv:2610.31989v1 [hep-th], 22 Oct 2026, prepared for JHEP. PDF: `papers/thong_2026.pdf`.

Equation numbers below are those of arXiv v1. Page numbers are the paper's printed numbers (PDF page = printed page + 2).

## In one paragraph

Integrate the centre position out of a CFT bilocal. The pair Casimir then acts only on the relative coordinate, and writing that coordinate as $x=e^\tau \mathbf n$ turns it into a cylinder Laplacian on $\mathbb R\times S^{d-1}$, the conformal boundary of AdS$_{d+1}$. A conformal family of dimension $h$ becomes a cylinder energy $\omega=h-d/2$. In a melonic CFT the reduced ladder kernel is exactly the Neumann boundary propagator $G_N$ of one free AdS scalar, so the geometric ladder series equals the resolvent that a quadratic boundary term produces. The ladder sum is therefore the Robin boundary propagator of that scalar: the constituent dimension $\Delta_\phi$ sets the mass, and the ladder coupling $g$ sets the Robin coefficient $\kappa$.

## Main question

Can one local free scalar in AdS$_{d+1}$ encode the conformal data of a melonic ladder, when the usual dictionary assigns one bulk field to each single-trace primary? (Sec. 1, p. 1)

## Physical system

Conformal melonic large-$N$ CFTs in $d$ dimensions, at leading order. Two examples:

| | SYK (Sec. 4.1) | Long-range Amit–Roginsky (Sec. 4.2) |
|---|---|---|
| Dimension | $d=1$ | $d=2$ |
| Disorder | yes | no |
| Constituent | Majorana fermion, $\Delta=1/q$ | real scalar, $\Delta_\phi=2/3$ |
| Bulk | AdS$_2$, odd sector ($J=1$) | AdS$_3$, even spins |
| $\nu$ | $1-2/q$ | $2/3$ |
| $m^2$ | $\nu^2-1/4$ (massless at $q=4$) | $-5/9$ |
| $\kappa$ | $g\,4^{-\nu}(q-1)$; $g=1$ is SYK (C.9) | $-3g/(8\pi)$ (4.5) |
| Stable range | $q>2$, $0\le g\le 1$ | $0\le g<g_c=1.02236136\ldots$ (4.7) |

## Important definitions

| Symbol | Meaning | Eq. |
|---|---|---|
| $X,\ x$ | pair centre $(x_1+x_2)/2$ and relative coordinate $x_{12}$ | Sec. 2 |
| $a$ | $2\Delta_\phi-d/2$, the weight in $W$ | Sec. 2 |
| $\nu$ | $d-2\Delta_\phi=d/2-a$ | (3.1) |
| $W$ | $(Wf)(\tau,\mathbf n)=e^{a\tau}f(e^\tau\mathbf n)$, removes the first derivative in $\tau$ | Sec. 2 |
| $\omega$ | $h-d/2$; the shadow $h\to d-h$ flips its sign | (2.6) |
| $F_g$ | ladder $(1-gK)^{-1}C_0$, with $K=C_0\mathcal V$ and $\mathcal V(y)=v_0|y|^{-2\nu}$ | (3.1) |
| $F_g^{(0)}$ | centre-integrated ladder | (3.2) |
| $c_\partial$ | seed normalisation | (A.6) |
| $\mathcal D_a$ | $[2\cosh(\tau-\tau')-2\mathbf n\cdot\mathbf n']^{-a}$ | (3.4) |
| $A,\ B$ | boundary coefficients: $\Phi\sim A\zeta^{\Delta_-}+B\zeta^{\Delta_+}$, with $\zeta=\cos\theta$, $\Delta_\pm=d/2\pm\nu$ | Sec. 3.2 |
| $\mathcal R_J(\omega)$ | response $B/A$ for the regular radial solution | (B.2) |
| $\kappa$ | Robin coefficient, $B=\kappa A$; identified as $\kappa=-gc_\partial/(2\nu)$ | Sec. 3.2 |

**Notation trap:** $J$ means the harmonic label in most of the paper but the SYK coupling in Sec. 4.1. At $d=1$, $J=0,1$ label even and odd functions on $S^0$, not spin. Also, $\mu(h)$ in (C.4) is unrelated to the $\mu$ of the N=2 draft.

## Main assumptions

- Leading large $N$, conformal limit. Higher-point functions and finite-$N$ effects are explicitly **not** identified (p. 2).
- A **single** ladder sum. Models with several ladders are outside the derivation; complex and SUSY SYK are named as cases where "we expect the bulk description to extend beyond the single scalar" (p. 6).
- Validity range (Sec. 5, p. 9; App. B): $0<\nu<1$, $a+J>0$ in every retained Neumann sector, and a positive Robin spatial operator. Outside it, "meromorphic continuation of the response does not by itself establish a stationary free-field interpretation."
- Centre integrals: absolutely convergent for $d/4<\Re\Delta_\phi<d/2$, otherwise defined by continuation. In SYK with $q\ge4$, the two Wick-exchange integrals diverge separately and must be combined before integrating (Sec. 3.1).
- The ladder is summed where the series converges, then continued as a resolvent (App. A).

## Main analytical results

The argument has five steps:

1. **Cylinder from the Casimir.** At zero centre momentum the pair Casimir reduces to $C_{\rm rel}=D_{\rm rel}(D_{\rm rel}-d)+\mathbf L^2$ (2.4). After $W$ it is $\partial_\tau^2+\mathbf L^2-d^2/4$ (2.5). The Casimir survives the centre integral but the full pair-generator action does not (A.2, A.4).
2. **Reduced kernel is a Neumann propagator.** $WK^{(0)}W^{-1}=c_\partial G_N$ with $G_N\propto P_\sigma\mathcal D_a$ (3.4). The mass $m^2=\nu^2-d^2/4$ (3.5) follows from this match. Geometry alone does not fix it.
3. **Robin scalar.** The boundary term $\int(\nu\kappa A^2-jA)$ gives $\widehat G_R=1/(2\nu[\kappa-\mathcal R_J(\omega)])$ (3.8), and so the double-trace resolvent $G_R=(G_N^{-1}+2\nu\kappa)^{-1}$ (3.9).
4. **Ladder = Robin.** With $\kappa=-gc_\partial/(2\nu)$:
   $$v_0(|x||y|)^aF_g^{(0)}(x,y)=c_\partial\,G_R(\tau-\tau',\mathbf n,\mathbf n';\kappa)\qquad\text{(3.11)}$$
   Equivalently, $k/(1-gk)=c_\partial\widehat G_R(h-d/2,J;\kappa)$ (3.12), with $k(h,J)=-c_\partial/(2\nu\mathcal R_J(\omega))$ (A.8). The spectrum is $\mathcal R_J(\omega)=\kappa$. Each rung is one $G_N$ joined at a boundary vertex (Fig. 3).
5. **Conformal data survive.** The same $c_m^2$ multiplies each family before and after reduction (3.13). The coefficients are unique if the reduced families are independent. This is verified by exponential asymptotics for SYK (C.2) and by explicitly inverting the centre integral for AR (D.3).

Further results:

- **Symmetry.** Nonzero $\nu\kappa$ breaks $2d$ AdS generators, because $A$ and $B$ have different weights (B.5). Only $\mathbb R\times SO(d)$ survives. Full AdS symmetry is not needed to encode the conformal data (p. 10).
- **Stability.** Stable for $\kappa>\mathcal R_J(0)$ in the lowest retained sector (App. B). In AR, $h_\pm=1\pm\omega_{00}$ collide at $g_c$ and become complex beyond it, though the mass stays above BF: this instability comes from the boundary condition.
- **SYK $h=2$.** Keep $g$ variable during the reduction and take $g\to1$ afterwards. This gives a finite mode at $\omega=3/2$, not a scalar zero mode (C.10). The divergent term drops out because $\mathcal I(2)=0$, i.e. $\mathcal P\Psi_2=0$ (C.12). Finite-coupling $h=2$ corrections are not included.
- **Relation to other constructions.** The DGJS interval problem (App. E) and the folded string (App. F) map to the same radial problem. These identify spectral problems, not spacetimes or interacting theories (p. 10).

## Important equations for the N=2 extension

$$\mathcal R_J(\omega)=\frac{\Gamma(-\nu)}{\Gamma(\nu)}\prod_\pm\frac{\Gamma[(J+\Delta_+\pm\omega)/2]}{\Gamma[(J+\Delta_-\pm\omega)/2]}\quad\text{(B.2)},\qquad \omega^N_{nJ}=a+J+2n\quad\text{(App. B)}$$

$$\mathcal J_h(\eta)\equiv\int_{-\infty}^\infty du\,\mathcal B_h\!\left(\frac{2}{1+\eta+2u^2}\right)=\frac{\pi}{2\omega\mathcal N(h)}e^{-\omega\,\mathrm{arcosh}\,\eta}\quad\text{(C.6)}$$

$$(\mathcal P\mathbb B_h)(\tau)=\mathcal I(h)\,e^{-(h-1/2)|\tau|},\quad \mathcal I(h)=\frac{1}{2\mu(h)\mathcal N(h)},\quad \Re h>\tfrac12\quad\text{(C.8)}$$

The condition $\Re h>1/2$ in (C.8) matters directly for the N=2 even sector: see [`docs/derivations.md`](../../docs/derivations.md), Sec. 3.

## Numerical methods

None in the usual sense: "No new data were generated or analysed" (p. 12). The results are analytic. Fig. 4 plots exact mode functions. The author acknowledges AI tools for figures and manuscript editing (p. 12).

## Relevant figures

- **Fig. 2** (p. 3): bilocal pair → relative space → cylinder.
- **Fig. 3** (p. 6): the rung-by-rung map between ladder and Robin series. This is the clearest statement of the idea.
- **Fig. 4** (p. 9): effective potentials and exact modes for SYK $q=4$, $q=2^{10}$ and AR $J=0,2$.

## Limitations

- A four-point statement only. It is explicitly not an equivalence of the full theories.
- Multiparticle scalar states are not multi-trace operators: energies $\omega=h-d/2$ do not add like dimensions (p. 9).
- "Radial quantisation" refers to the reduced pair Casimir, not a full state–operator map (p. 9).
- SYK's soft-mode physics is not captured: the Schwarzian match and a fluctuating boundary are left open (p. 11).
- The bulk scalar's spacetime has no dynamical gravity, and the $\tau$ of the cylinder is $\log|T_{12}|$, not microscopic time.
- Convergence relies on continuation prescriptions, which must be used consistently across compositions (App. A).

## Independent checks

Run with `scripts/check_robin_identities.py`; record in `results/data/robin_identity_checks/`. These test internal consistency of the printed formulas, not the whole derivation.

| Statement | Result |
|---|---|
| $k_c(h)\,\mathcal R_-(h-\tfrac12)=4^{-\nu}(q-1)$ (App. C.1), $q=4,6,10$, real and complex $h$ | agrees (see caveat) |
| $c_\partial=-2\nu4^{-\nu}(q-1)$ from (C.3), (A.6) | agrees |
| $q=4$ spectrum from $\mathcal R_-(\omega)=\kappa$: $h=2,\ 3.77354,\ 5.67946,\ 7.63197$ | agrees with roots of $k_c=1$ |
| $G_N$ is the Neumann response, (3.4): $d=1$ odd ($\nu=0.3,0.8$) and $d=2$, $J=0$ (AR), by direct Fourier integral | agrees |
| Block centre integral (C.6) | agrees |
| Coupling identity (E.3) | agrees |
| AR $g_c=1.02236136077$ (4.7); $c_\partial=1/(2\pi)$ (A.5, A.6, D.1) | agrees |

By hand: (2.5), (2.6), the second equality of (A.8), the sign pattern and beta functions of (C.3), and the trigonometric reduction in (C.8) all check out. The AdS$_3$ wave operator in (4.6) also follows from its metric.

**Caveat:** the $k_c$ used here was transcribed from memory of Maldacena–Stanford (3.35), which is not in `papers/`. It passes $k_c(2)=1$ and the identity above, but should be compared to the source.

**Not checked:** the Wick combinatorics of $v_0$; (C.10)–(C.12) beyond $\mathcal I(2)=0$; the absolute-convergence and term-by-term arguments of App. C.2 and D.1; the OPE residues; and the eigenvalue factor of 1/2 relative to Benedetti–Delporte's printed (6.19) (App. D).

## Relationship to our project

This is the construction we are extending. Its stated limits show where the work is:

- It handles a **single** ladder with **one** constituent dimension. N=2 SYK has a mixed fermion–boson neutral channel and a separate charged channel.
- Its validity range is **not automatically met** in N=2. The auxiliary boson has $\nu_b=-\delta<0$, and the even-parity fermion channel has $a+J=-\mu<0$.
- Its first step (cylinder from the Casimir) is the part with no N=2 version yet. The draft starts from the reduced kernels, not from the $su(1,1|1)$ Casimir.
- Its SYK centre-reduction analysis (App. C.2) has no N=2 counterpart yet, and that is what would justify carrying OPE data through the reduction.

Details are in [`docs/derivations.md`](../../docs/derivations.md) and [`docs/research_questions.md`](../../docs/research_questions.md).

## Update after the literature review (2026-11-06)

- **Cited Maldacena–Stanford numbers** are from a later version than the v1 in `papers/`. For example, Thong's [MS](3.35) is v1 (3.73). The kernel itself agrees.
- **(C.8) has an even-sector analogue,** with $\tan(\pi h/2)\to\cot(\pi h/2)$, found by quadrature. See `docs/derivations.md` 7.1.
- **(C.8) has a pole at $h=1$** that Thong never meets, as the physical roots there stay in $(2n+1+2\Delta,2n+2)$. N=2 does meet it. The family there is a step function on $\chi>1$, `docs/derivations.md` 7.2.
- **His p. 6 remark on complex SYK** can be sharpened at zero charge: one scalar with a parity-dependent $\kappa$, `docs/derivations.md` Sec. 9.

