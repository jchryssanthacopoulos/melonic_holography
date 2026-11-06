# Fu, Gaiotto, Maldacena, Sachdev 2016

**Citation.** W. Fu, D. Gaiotto, J. Maldacena, S. Sachdev, *Supersymmetric SYK models*, Phys. Rev. D 95 (2017) 026009, arXiv:1610.08917. PDF: `papers/fu_gaiotto_maldacena_sachdev_2016.pdf`. Tag: **[FGMS]**.

**Read:** Secs. V (N=2), VI (four-point function and spectrum), VII (conclusions). Secs. I–IV (N=1) were only searched, and need a proper read before the N=1 stepping stone.

## Main question

Define supersymmetric SYK, find its low-energy structure and spectrum, and identify the soft modes.

## Physical system

N=2: complex fermions with $Q=i\sum C_{ijk}\psi^i\psi^j\psi^k$, $H=\{Q,\bar Q\}$ (5.1)–(5.3), $\overline{|C|^2}=2J/N^2$ (5.4). Generalizes to $\hat q$, with up to $2\hat q-2$ fermions in $H$. Note that (5.4) has a different normalization of $J$ from [PSV](2.2).

## Main assumptions

Large $N$, IR superconformal solution, and a soft-mode effective action fixed by symmetry.

## Main analytical results

- **Dimensions and charges**: $\Delta_\psi=1/(2\hat q)$, R-charge $1/\hat q$, so $\psi$ is a chiral primary with $Q_R=2\Delta$. Global symmetry $\mathbb Z_{\hat q}$.
- **Supersymmetry is unbroken**: the Witten index (5.5), (5.6) saturates the large-$N$ ground state entropy.
- **Superspace** (5.8)–(5.11), including the chiral field and $\Delta_{12}$. **Super-reparameterizations** (5.14) and $SU(1,1|1)$ (5.23)–(5.26).
- **Super-Schwarzian** (5.27)–(5.29). Bosonic part (5.37): Schwarzian plus $-2(\partial_\tau\sigma)^2$, with the relative coefficient fixed by N=2, so low-energy thermodynamics has one coefficient $\alpha_s$ (5.38)–(5.40). Density of states (5.42), winding sum (5.43).
- **Spectrum** (6.2), (6.3). Two towers with $h=\tfrac32,3.3211,5.2409,\ldots$ and $h=\tfrac32,3.5659,5.5949,\ldots$, each with $h\pm\tfrac12$ partners.
- **Extra IR symmetry** (6.4)–(6.6): $G_{\psi\bar\psi}\to[f'(\tau_1)h'(\tau_2)]^\Delta G(f,h)$ with $f$ and $h$ independent. It is broken by the UV boundary condition unless $h=f$.

## The statement that matters most for us

Sec. VI, p. 31: "we expect that by supersymmetry ... the multiplet coming from the symmetric tower with dimensions $(1,3/2,2)$ should lead to the superschwarzian while the second one, coming from $h_{a,m}$, also with dimensions $(1,3/2,2)$ should correspond to operators in the IR theory. As before these arise from symmetries of the low energy equations, namely (5.21)."

So **two different $H=1$ multiplets**, with different physics:

| Multiplet | [FGMS] | Project | Physics |
|---|---|---|---|
| "$s$" | $k_s$ tower | $S$-branch | Enhanced; divergent in the conformal ladder; super-Schwarzian |
| "$a$" | $k_a$ tower | $A$-branch | Not divergent; "operators in the IR theory"; from $f\neq h$ |

The label map was checked: [FGMS](6.2) gives $k_s(h)=k_S(h-\tfrac12)$, $k_a(h)=k_A(h-\tfrac12)$.

They state this as an expectation: "we have not filled out all the details".

## Numerical methods

Exact diagonalization (up to $N=28$ for N=1, and degeneracies for N=2 at $\hat q=3$, (5.7)) against the large-$N$ solution.

## Limitations

The four-point function is not computed. $\hat q=3$ is used for most formulas. The zero-mode determinant is not exact in N=2 (Sec. V end).

## Relationship to our project

- **Spectrum**: [N2](15) agrees with (6.2), numerically for $q=3,5,7$.
- **C2** (the $A$-branch $H=1$ modes): the quoted statement supplies the physics. Together with the result in `docs/derivations.md` Sec. 7, the tension is resolved.
- **C5** (the physical $H=1$ coefficient): the super-Schwarzian (5.27)–(5.38) is what's needed to go beyond the draft's finite part [N2](33). Not yet done.
- **$J$ conventions differ**: (5.4), [PSV](2.2) and [B18](3.2) must not be mixed.
