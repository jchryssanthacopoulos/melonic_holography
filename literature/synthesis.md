# Literature synthesis

The project's current understanding of the literature, updated as papers are reviewed. Per-paper detail is in [`notes/`](notes/); references are in [`bibliography.md`](bibliography.md). Eleven papers are in `papers/`. **Several were read only in part**, with the coverage stated in each note and in the bibliography. Statements about any other paper are second-hand and are marked.

*Last updated 2026-11-06.*

## The problem

Melonic large-$N$ CFTs are solvable: propagators resum iterated melons, and leading connected four-point functions are ladder sums. In the usual large-$N$ dictionary each single-trace primary corresponds to one bulk field, so the infinite bilinear tower of SYK seems to require an infinite tower of fields. Thong asks if one local free field can encode the same data.

## The solvable structure, from SYK to N=2

The same recipe runs through the literature, from [Maldacena–Stanford](notes/maldacena_stanford_2016.md) to [complex SYK](notes/bulycheva_2017.md), [N=1](notes/murugan_stanford_witten_2017.md), and [N=2 in components](notes/peng_spradlin_volovich_2017.md) and [in superspace](notes/bulycheva_2018.md): find the two-point function; diagonalize the kernel on three-point structures; expand the seed in Casimir eigenfunctions; sum a geometric series; close a contour.

| | $\mathcal N=0$ | complex | $\mathcal N=1$ | $\mathcal N=2$ |
|---|---|---|---|---|
| Casimir eigenvalue | $h(h-1)$ | $h(h-1)$ | $h(h-\tfrac12)$ | $h^2$ |
| Shadow | $1-h$ | $1-h$ | $\tfrac12-h$ | $-h$ |
| Inner product | positive | positive | no Hilbert space | $d\chi/[\chi(1-\chi)]$, not positive |
| Seed in the space? | yes | antisym. yes, sym. no | | marginally, by p.v. |
| Protected modes | $h=2$ | $h=2$; $h=1$ | $h=\tfrac32$ multiplet | $(1,\tfrac32,\tfrac32,2)$, twice |

N=2 is special in three ways. Its kernel needs both symmetric and antisymmetric sectors. It isn't even in $h$: $k^A(h)=k^S(-h)$. And chirality makes the four-point function purely bosonic in $\chi$.

**Two $H=1$ multiplets.** [Fu et al.](notes/fu_gaiotto_maldacena_sachdev_2016.md) expect one to give the super-Schwarzian, and the other IR operators tied to an extra symmetry. [PSV](notes/peng_spradlin_volovich_2017.md) confirm that only one diverges. The super-Schwarzian contribution to the four-point function has not been computed in any paper we have read.

**A line of theories.** [Gross–Rosenhaus](notes/gross_rosenhaus_2017.md) replace the kinetic term by a bilocal one, giving a line of CFTs with $k\to(1-2b)k$, a free end and an SYK end. Their bulk is a tower of fields on rigid AdS$_2$.

## Thong's construction

Thong's central claim (arXiv:2610.31989, Eq. 3.11): after integrating out the pair centre, the melonic ladder sum **is** the Robin boundary-to-boundary propagator of a single free scalar in rigid AdS$_{d+1}$, whose boundary is the relative cylinder.

| CFT side | Bulk side |
|---|---|
| Constituent dimension $\Delta_\phi$ | Mass, $m^2=\nu^2-d^2/4$ with $\nu=d-2\Delta_\phi$ |
| Rung kernel | Neumann boundary propagator $G_N$ |
| Ladder coupling $g$ | Robin coefficient $\kappa=-gc_\partial/(2\nu)$ |
| Ladder sum | Robin propagator, a double-trace resolvent |
| Family dimension $h$ | Cylinder energy $\omega=h-d/2$ |
| Kernel equation $k=1/g$ | $\mathcal R_J(\omega)=\kappa$ |
| Complex dimensions in AR | Boundary-condition instability above BF |

Three aspects stand out:

1. **The coupling is in the boundary condition**, not in a potential.
2. **AdS symmetry is broken and it doesn't matter.** Robin breaks $2d$ generators, leaving $\mathbb R\times SO(d)$, yet dimensions and OPE coefficients are retained.
3. **The earlier constructions are related.** Per Thong, the DGJS and folded-string problems map to the same radial problem, as identifications of spectral problems, not of spacetimes. We haven't read those papers.

It is deliberately modest: four-point data only. Its foundation for stability is [Ishibashi–Wald](notes/ishibashi_wald_2004.md), which classifies boundary conditions for one field but says nothing about matrix conditions. The closed-form identities all hold (see [`notes/thong_2026.md`](notes/thong_2026.md)). Its analytical arguments for convergence and term-by-term reduction haven't been independently checked.

Thong's only remark on N=2 is that for "richer constituent field content and multiple ladder sums, such as complex or SUSY SYK ... we expect the bulk description to extend beyond the single scalar" (p. 6).

## Supersymmetric boundary conditions

[Correa et al.](notes/correa_giraldo-rivera_silva_2019.md) classify supersymmetric mixed conditions for massless fields on AdS$_2$. A plain Robin condition keeps 2 of 4 supercharges and breaks scale invariance in Poincaré time. A derivative condition on the scalar with an algebraic one on the fermion keeps all 4. The AdS$_2$ N=2 multiplet needs $m_B^2=m_F^2-m_F$, and the Killing spinors have global energy $\pm\tfrac12$. **Nothing here covers massive fields, which is the case the draft needs.**

## The project's N=2 draft

`research/pdf/n2_susy_generalization_note.pdf` is AI-generated. It proposes two scalars with a constant *matrix* Robin condition for neutral channels; a Dirac response for charged channels; and a fixed Abel-type transform that rebuilds the superconformal four-point function. It claims a four-point representation only.

**Assessment after comparing with the sources** (details in [`../docs/derivations.md`](../docs/derivations.md) Sec. 5):

- **Its physical input is right.** Kernels, spectra, blocks, seed and densities agree with [PSV], [Bulycheva 2018] and [Fu et al.].
- **Its claimed misprints are real.** The two sign problems in [Bulycheva 2018] are confirmed, and we found further ones, including a factor of 2 in the printed four-point function.
- **Its central identity holds where tested.** Computing the centre-integrated four-point function directly from [PSV] with Thong's (C.8) gives exactly the auxiliary covariance, at finite $g$, in both parities, in the fermion entry, for $\omega<4$.

## New in this project

These are project-derived, not from the literature, and **not externally reviewed or searched for in prior work** (D5).

- **The relative cylinder from the super-Casimir.** At zero total supermomentum $C\to L_0^2$, so $\omega=H$ and shadow $-H$. The relative space is bosonic. In general $d_s=1-\mathcal N/2$.
- **Why the draft has a matrix.** Ordinary centre integration keeps both bottom and top components of a multiplet. Super-centre integration would keep one, so a single Dirac field is the conjectured N=2 analogue of Thong's single scalar.
- **Exceptional families.** The $A$-branch $H=1$ modes are finite step functions on $\chi>1$, invisible to the $0<\chi<1$ OPE and to the draft's reconstruction. Their centre integral gives exactly the mode, via the Legendre generating function. This resolves an apparent contradiction with Thong's (3.13).
- **A supersymmetric line.** A superspace bilocal kinetic term gives a line with uniform $k\to gk$, on which the draft's formal $g$ and threshold $g_*$ may be physical.
- **Complex SYK is one scalar.** Its symmetric channel is exactly Thong's $J=0$ sector, with a different $\kappa$.
- **Stability.** The matrix boundary condition is self-adjoint with $\omega^2\ge0$ for all real $g$. In the even sector, the continued seed differs from the stationary covariance by one homogeneous mode, with a formal zero mode at $g_*<1$.

## What's still missing

The super-centre integral of a superblock, which would make the Dirac statement a theorem. The surviving symmetry. The super-Schwarzian coefficient. Any comparison with finite-$N$ data. A literature search.
