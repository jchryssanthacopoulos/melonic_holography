# Derivations and theoretical notes

This is the project's working record of what is derived and how far it is verified.

**Status labels**

- **Paper claim**: stated in a primary paper, not independently verified unless noted.
- **Draft claim**: stated in `research/pdf/n2_susy_generalization_note.pdf` (AI-generated, unreviewed).
- **Checked**: the formula was tested by hand or by a script in `scripts/` (records: `results/data/robin_identity_checks/`, `results/data/literature_checks/`). Unless it is compared against a source paper, this shows the formula is right as mathematics, not that it describes N=2 SYK.
- **Derived here**: derived in this project, with the argument shown. Not externally reviewed.
- **Open**: not established.

Equation tags: **[T]** is Thong 2026 arXiv v1; **[N2]** is the draft; **[PSV]**, **[B18]**, **[B17]**, **[FGMS]**, **[MS]** (v1 numbering), **[GR]**, **[CGS]**, **[IW]**, **[MSW]**, **[GKST]** are as in `literature/notes/`.

**Where things stand (2026-11-06).** The draft's physical input agrees with its sources (Sec. 5). B1 now has a derivation (Sec. 6). The $A$-branch tension C2 is resolved (Sec. 7). D2 and D1 have new answers (Secs. 8, 9). The largest remaining gaps are in Sec. 4.

## 1. Notation

Follow Thong's symbols wherever they apply. The draft's additions, and the name collisions to watch for:

| Project symbol | Meaning | Thong | Draft |
|---|---|---|---|
| $\delta$ | $1/q$, with $q$ odd the supercharge order | (absent) | $\delta$ |
| $\Delta_\psi,\ \Delta_b$ | $\delta/2,\ (1+\delta)/2$ | $\Delta_\phi$ | $\Delta_\psi,\ \Delta_b$ |
| $\nu$ | fermion index $1-\delta$ | $d-2\Delta_\phi$ with $\Delta_\phi=\Delta_\psi$ | $\nu$ |
| $\mu$ | $\tfrac12-\delta$; the Dirac mass | **not** $\mu(h)$ of [T](C.4) | $\mu$ |
| $a_\psi,\ a_b$ | $-\mu,\ 1-\mu$ | $a$ | $a_\psi,\ a_b$ |
| $\tau$ | cylinder time | $\tau$ | $s$ (cylinder), $t$ (reconstruction) |
| $\sigma,\ p$ | $\sigma=(-1)^p$; $p$ is the fermion-pair parity | $J$ at $d=1$ labels parity, and $\sigma$ the retained pair parity | $\sigma,\ p$ |
| $h,\ H$ | component dimension, superprimary dimension | $h$ | $h,\ H$ |
| $\omega$ | $h-\tfrac12$ | $\omega$ | $\omega$ |
| $J$ | N=2 coupling (draft); harmonic label (Thong) | **both** meanings | coupling |

The bulk field $\Phi_b$ has nothing to do with Thong's $b$ (the SYK propagator amplitude), and $b_\psi$ is an amplitude.

**Labels across papers. These collide, so always check.** Kernel entries were verified in `scripts/check_literature_inputs.py`.

| Project | [PSV] kernels | [PSV] four-point | [B18] | [FGMS] |
|---|---|---|---|---|
| $p=1$ (odd fermion pair) | $k^{a,\pm}$ | $f^s$, $\Phi^s$ | | |
| $p=0$ (even fermion pair) | $k^{s,\pm}$ | $f^a$, $\Phi^a$ | | |
| $k_A(H)$ | $k^{a,+}(H)=k^{a,d}(H+\tfrac12)=k^{s,-}(H+1)$ | | $-k^A(H)$ as printed | $k_a(H+\tfrac12)$ |
| $k_S(H)$ | $k^{s,+}(H)=k^{s,d}(H+\tfrac12)=k^{a,-}(H+1)$ | | $-k^S(H)$ as printed | $k_s(H+\tfrac12)$ |
| $z$ | $\chi/(\chi-1)$, with $\chi<0$ | | $\chi$ | |
| $\delta=2\Delta_\psi$ | $1/q$ | | $2\Delta$ | $1/\hat q$ |

$J$ is normalized differently in [PSV](2.2) and [FGMS](5.4), and $b$ differs in [PSV](2.10) and [B18](3.9) (see Sec. 5.2).

## 2. The N=2 draft: structure and claim status

### 2.1 The construction

The draft works at leading large $N$, at the charge-symmetric saddle, in the conformal limit. The steps:

1. **Neutral channels mix.** $\psi\bar\psi$ (parity $p$) and $b\bar b$ (parity $1-p$) form a $2\times2$ problem, $\mathsf K=\mathsf C_0\mathsf V$ [N2](3). The powers cancel so that the boundary mixing is constant: $1-a_i-a_j-v_{ij}=0$.
2. **Two scalars.** $\Phi_\psi$ has $m^2=\mu(\mu+1)$ with alternate quantization. $\Phi_b$ has $m^2=\mu(\mu-1)$ with standard quantization. $X=(A_\psi,B_b)^T$ [N2](5), (7).
3. **Matrix Robin.** $\widehat{\mathsf k}_p=\mathsf D_p\mathsf M_p$ [N2](9), so $\mathsf G_{g,p}=(\mathsf D_p^{-1}-g\mathsf M_p)^{-1}$ [N2](10). The boundary conditions are [N2](11):
   $$-2\nu B_\psi=g\sigma L_q(2\mu A_\psi+2\delta B_b),\qquad A_b=g\sigma L_qA_\psi.$$
4. **Supersymmetric branches.** The eigenvalues are $\sigma L_qD_\psi(\mu\pm\omega)$, with frequency-dependent eigenvectors, so a constant rotation cannot decouple them. They reorganize into $k_A$ and $k_S$ [N2](15), with multiplets $(H,H+\tfrac12,H+\tfrac12,H+1)$ [N2](16), (17).
5. **Dirac field.** The charged channels use a Dirac field of mass $\mu$ with $R_{F,p}=2\nu R_{\nu,p}(\Omega-\tfrac12)/(\Omega-\nu)$ [N2](21), and $k_\alpha=\sigma L_q/R_{F,p}$ [N2](22). Its first-order radial equations [N2](19) tie the two scalars together.
6. **Reconstruction.** A fixed Abel-type transform $\mathcal T_2$ sends $e^{-Ht}$ to the superblock $\varphi_H$ [N2](25). A spectral density and contour rebuild the four-point function [N2](26), (27).
7. **Exceptional terms.** The $H=0$ residue [N2](30), and a finite-part prescription for the $H=1$ multiplet $(1,\tfrac32,\tfrac32,2)$ [N2](31)–(33).
8. **Negative result.** The fixed boundary domain does not preserve the AdS isometry algebra, so it cannot preserve all four AdS N=2 supercharges [N2](35).

### 2.2 What has been checked

All entries below pass, for $q=3,5,7$ where $q$ applies. Where the entry is numerical, residuals are at most $10^{-8}$, and typically far smaller.

| [N2] | Statement | How |
|---|---|---|
| (5) | masses from $\nu^2-\tfrac14$, $\delta^2-\tfrac14$; boundary powers | numeric, hand |
| (12) | $r_\psi$ equals the chain integral; $r_\psi/c_\psi=\ell_\psi$ | numeric |
| (13) | $R_{\nu,p}R_{\delta,1-p}=(\mu^2-\omega^2)/(4\delta\nu)$ | numeric, hand |
| (9), (14), (16) | eigenvalues of $\mathsf D_p\mathsf M_p$ equal $\{k_A,k_S\}$ with the stated shifts | numeric (trace and det, real and complex $h$) |
| (15) | $k_A(1)=k_S(1)=1$, $k_A(0)=-(q-1)$ | numeric, hand |
| (19) | squaring gives the two masses | hand |
| (20) | boundary maps from (19) | hand |
| (11) from (20), (23) | the Dirac condition implies the matrix condition | hand |
| (21), (22) and the top-scalar form | | numeric |
| Sec. 3 | $H_S=2.8211441710$, $H_A=3.0659046429$ at $q=3$ | numeric |
| (24), (25) | contiguous identity; $\mathcal T_2e^{-Ht}=\varphi_H$; the underlying Abel relation | numeric (one case re-confirmed in plain Python) |
| (29) | seed identity | numeric |
| (27) | at order $g^0$ gives $z^\delta$; at order $g^1$ gives the table of (34) | numeric (contour + integer sum) |
| (34) | direct convolution reproduces the printed table | numeric |
| (30) | $H=0$ residue | hand |
| (31), (33) | $b_1$, $b_2$, $a_0$, and the two $q=3$ coefficients | numeric |
| (35) | the transformed boundary mismatch | hand |
| Sec. 5.1 | large-$|{\rm Im}\,H|$ decay $O(|{\rm Im}\,H|^{4\delta-5/2})$ | hand, leading order |

### 2.3 What has not been checked

The checks in 2.2 show internal consistency. Comparison with the sources is now done in Sec. 5. What is left:

1. **Finite $g$ in the contour form.** (27) has only been tested as a contour through order $g^1$ here. Its residues (28) are now checked against [PSV] at finite $g$ (Sec. 5.4), which leaves the arcs at infinity, the explicit integer sum, and order $g^2$ against a direct two-rung integral.
2. **The entries of (12) other than 11.** Sec. 5.4 tests $\psi\bar\psi\psi\bar\psi$ only.
3. **The $H=1$ finite part** (32), (33) beyond arithmetic. The draft says a different subtraction changes the $\varphi_1$ coefficient, and that the physical value requires the super-Schwarzian.
4. **Term-by-term reduction** of the infinite family sum. Tests stop at $\omega<4$.
5. **The calculation records** `n2_susy_robin_extension.md` and `n2_completion_supplement.md` are cited by the draft but not in this repository.

### 2.4 Mismatch with the published paper

The draft was written against a pre-arXiv version of Thong's paper ("supplied draft", titled *Holography with a single bulk field in AdS$_{d+1}$*). Consequently:

- Its citations ("Eqs. (3), (11)–(13)") do not match the numbering of v1, and we don't have that version to map them.
- The "ordinary Abel map" $\mathcal A$ of [N2](25) **does not appear in v1**. The nearest statement is the forward integral [T](C.6). $\mathcal A e^{-(H-1/2)t}=F^{\rm bos}_H/2$ has the form of an Abel inversion of (C.6) with $\eta=2/z-1$. **Checked** numerically, but it should now be derived in this project rather than inherited.

## 3. The auxiliary problem against Thong's validity range

This section addresses one of the draft's disclaimers, which says it does not supply "a positive auxiliary vacuum in every continued sector" ([N2], Sec. 7).

**Summary**

| Thong's condition | N=2 fermion rail | N=2 boson rail |
|---|---|---|
| $0<\nu<1$ | yes: $\nu=1-\delta$ | **no**: $d-2\Delta_b=-\delta$. Handled with $\rho=\delta$ and standard quantization. |
| $a+J>0$ | $p=1$: $1-\mu>0$. **$p=0$: $-\mu<0$, fails.** | yes: $1-\mu+(1-p)>0$ |
| Positive Robin operator | see below | see below |

### 3.1 The matrix boundary condition is self-adjoint (Derived here)

**Assumptions.** Real $g$; both scalars have index in $(0,1)$ and so are limit-circle at the boundary; fixed parity on the strip, so that one endpoint carries independent data.

**Variables.** $c\equiv g\sigma L_q$, and boundary data as in [N2](5).

**Derivation.** For $\Phi=A\zeta^{1/2-\rho}+B\zeta^{1/2+\rho}$, the Wronskian at the boundary is
$$\Phi_1\partial_\zeta\Phi_2-\Phi_2\partial_\zeta\Phi_1\to2\rho(A_1B_2-B_1A_2).$$
The two-channel boundary form is therefore
$$\Omega=2\nu(A_{\psi1}B_{\psi2}-B_{\psi1}A_{\psi2})+2\delta(A_{b1}B_{b2}-B_{b1}A_{b2}).$$
Impose (11): $2\nu B_\psi=-c(2\mu A_\psi+2\delta B_b)$ and $A_b=cA_\psi$, leaving $(A_\psi,B_b)$ free. Then
$$\Omega=-cA_{\psi1}(2\mu A_{\psi2}+2\delta B_{b2})+c(2\mu A_{\psi1}+2\delta B_{b1})A_{\psi2}+2\delta c(A_{\psi1}B_{b2}-B_{b1}A_{\psi2})=0.$$
The boundary condition is Lagrangian, so $\mathrm{diag}(-\partial_\theta^2+m_\psi^2\sec^2\theta,\ -\partial_\theta^2+m_b^2\sec^2\theta)$ is self-adjoint on it in the ordinary positive $L^2$ norm. **So $\omega^2$ is real for all real $g$**, even though the matrix eigenvectors depend on frequency.

**Check.** $g=0$ gives $B_\psi=A_b=0$, the two decoupled free problems.

### 3.2 There are no tachyons (Derived here)

**Claim.** For real $g\neq0$, no mode has $\omega^2<0$.

**Derivation.** Put $\omega=iy$ with real $y\neq0$. $R_{\nu,p}(iy)$ is real, finite and nonzero: its gamma factors come in conjugate pairs with nonzero imaginary part. By [N2](13), $R_{\delta,1-p}(iy)$ is too, so $\mathsf D_p$ is invertible. A mode requires an eigenvalue of $\widehat{\mathsf k}_p$ to equal the real number $1/g$. The eigenvalues are $\sigma L_qD_\psi(\mu\pm iy)$ with real nonzero $D_\psi$, and their imaginary part is nonzero. Contradiction. $\square$

Combined with 3.1, the spectrum satisfies $\omega^2\ge0$. **This differs from AR**, where $\omega^2$ turns negative beyond $g_c$ ([T] Sec. 4.2).

### 3.3 The even sector has a formal zero mode at $g_*<1$ (Derived here + numerics)

At $\omega=0$ both eigenvalues of $\widehat{\mathsf k}_p$ coincide, with $\widehat{\mathsf k}_p$ non-diagonalizable. The zero-mode condition is $1/g=k_S(\tfrac12)$ for $p=0$ and $1/g=k_A(\tfrac12)$ for $p=1$.

| $q$ | $k_S(\tfrac12)$ | $g_*=1/k_S(\tfrac12)$ | $k_A(\tfrac12)$ |
|---|---|---|---|
| 3 | 4.59057 | 0.21783783 | −0.329588 |
| 5 | 2.87947 | 0.34728561 | −0.747558 |
| 7 | 2.53261 | 0.39484963 | −0.999908 |

Because $k_A(\tfrac12)<0$ in these cases, $p=1$ has no crossing at $g>0$. In $p=0$ the lowest mode starts at $|\omega|=\mu$ at $g=0$, falls to $0$ at $g_*$, and rises to $\omega=\tfrac12$ at $g=1$. Since $\omega$ passes linearly through zero, $\omega^2$ **touches zero and returns** instead of going negative. That is consistent with 3.2.

Near $g_*$ the numerics show the $2\times2$ residue of $\mathsf G_{g,0}$ growing like $1/\omega$, as the frequency-dependent projectors in [N2](18) would suggest (e.g. $q=3$: $\omega\times$ weight $\approx0.26$ at $g=0.25$, $0.29$ at $g=0.5$, $0.30$ at $g=0.9$). So $g_*$ is a threshold in the sense of [T] App. B, with no stationary vacuum exactly at that point.

$g$ is formal in N=2 and the physical point is $g=1$, so this threshold is below where the physics lives, as long as nothing requires continuity in $g$ from $g=0$. **See 3.5.**

### 3.4 The continued seed is not the stationary covariance when $a+J<0$ (Derived here + numerics)

For $p=0$ the reduced fermion rail grows at large $|\tau|$, $\mathcal D_a\sim e^{\mu|\tau|}$, so it is not a stationary covariance. The precise relation is as follows. **Checked** numerically for $q=3,5$ at two values of $\tau$: the reduced kernel of [T](3.4) equals the mode sum
$$\frac{\Gamma(a)}{2\sqrt\pi\,\Gamma(1-\nu)}\big(\mathcal D_a^{\rm same}+\sigma\mathcal D_a^{\rm opp}\big)=\sum_{n\ge0}w(\omega_n)e^{-\omega_n|\tau|},\quad\omega_n=a+J+2n,\quad w(\omega)=\frac{1}{2\nu\mathcal R_J'(\omega)},$$
where $\mathcal D_a^{\rm same},\mathcal D_a^{\rm opp}$ are $\mathcal D_a$ at $\mathbf n\cdot\mathbf n'=\pm1$. The $n=0$ term keeps the continued exponent. When $a+J=-\mu<0$, it has $w(-\mu)=-w(\mu)<0$, so it is a *negative-coefficient growing* exponential. The stationary Neumann covariance of the same scalar has the *positive decaying* mode $w(\mu)e^{-\mu|\tau|}$ instead. These have the same meromorphic Fourier response, and
$$G^{\rm continued}-G^{\rm stationary}=-2w(\mu)\cosh(\mu\tau),$$
which solves the homogeneous equation. So the seed fails Thong's condition by exactly one homogeneous mode. $q=3$: $w(\mu)=1.42635$. $q=5$: $w(\mu)=0.531746$.

This matches [T](C.8), which requires $\Re h>\tfrac12$ for the centre integral to converge. The free $p=0$ pole has $h=\delta<\tfrac12$. The draft's remark that "the separate fermion rail has an infrared divergence" is the same fact.

### 3.5 Status at $g=1$ (numerics; Open in part)

At $g=1$, the lowest roots of $\det(1-g\widehat{\mathsf k}_p)$ for $\omega<4$:

| $q$ | $p$ | $\omega$ |
|---|---|---|
| 3 | 0 | 0.5, 1.5, 2.32114, 3.5659 |
| 3 | 1 | 0.5, 1.5, 2.5659, 3.32114 |
| 5 | 0 | 0.5, 1.5, 2.05953, 3.56561 |
| 5 | 1 | 0.5, 1.5, 2.56561, 3.05953 |
| 7 | 0 | 0.5, 1.5, 1.91414, 3.5542 |
| 7 | 1 | 0.5, 1.5, 2.5542, 2.91414 |

At $q=3$, the two higher values in each row are $H\mp\tfrac12$ for the printed $H_S,H_A$ of [N2]; the draft prints no roots for $q=5,7$. All have $\omega>0$. At each root listed for $q=3,5,7$ and $g\in\{0.25,0.5,0.9,1\}$, $-\mathrm{Res}\,\mathsf G_{g,p}$ is **positive semi-definite of rank one**, within rounding. That is what a stationary free-field spectral sum needs.

**Conclusion (provisional).** At physical $g=1$ the auxiliary system looks like a healthy self-adjoint problem, with every mode above threshold. That is stronger than the draft claimed. The limitations:

- The root scan stops at $\omega<4$ and tests finitely many $(q,g)$.
- **Mostly settled, see Sec. 5.4:** the reduced $g=1$ ladder, computed directly from [PSV]'s four-point function with [T](C.8), has exactly the stationary covariance's modes and weights, at least in the 11 entry and for $\omega<4$. No continuity in $g$ through $g_*$ is needed, since [PSV]'s contour at $\Re h=\tfrac12$ defines the $g=1$ answer directly.

### 3.6 A possible tension: A-branch modes at $H=1$ (Resolved; see Sec. 7)

*Historical. Kept as first stated.*

At $g=1$ both branches have $H=1$ roots, since $k_A(1)=k_S(1)=1$:

| Branch | Sector, mode | In the draft's four-point function |
|---|---|---|
| $S$ | $p=0,\ \omega=\tfrac12$ ($h=1$); $p=1,\ \omega=\tfrac32$ ($h=2$) | the $(1,\tfrac32,\tfrac32,2)$ multiplet |
| $A$ | $p=1,\ \omega=\tfrac12$ ($h=1$); $p=0,\ \omega=\tfrac32$ ($h=2$) | weight vanishes: $c_A\propto\cot(\pi H/2)\to0$ [N2](28) |

Numerically, all four modes have nonzero positive residue in $\mathsf G_{1,p}$ (e.g. $q=3$: 0.592, 0.493, 0.246, 0.205). If [N2](12) holds with these weights, the reduced component covariance contains the $A$-branch modes, while the chiral–antichiral reconstruction drops them. Thong's (3.13) says the same $c_m^2$ multiplies each family before and after reduction. So this needs one of: these modes appear in other orderings or components; or (12) and (27) use different data. The draft's phrase "in this particular seed and ordering" allows the first. This has to be resolved.

## 4. Gaps relative to Thong's derivation

Thong's construction has steps the draft takes as given.

1. **Cylinder from a Casimir. Now derived,** Sec. 6: the eigenvalue is $H^2$ because $d_s=0$, with no odd relative coordinate.
2. **Centre reduction of families. Done for ordinary centre integrals of the component blocks,** in both parities, Sec. 7. **Open:** the super-centre integral of a superblock, for which the draft's $\mathcal T_2$ would presumably be the inverse. The draft calls $\mathcal T_2$ "a prescribed reconstruction, not an inverse on arbitrary center-integrated data".
3. **Independence and term-by-term reduction** of families ([T] App. A, C.2). Open. The new step-function families (Sec. 7) need including.
4. **Convergence.** The fermion rail is IR divergent and the boson rail UV divergent; the draft uses "a common homogeneous-distribution continuation". Thong requires the prescription to allow centre shifts and interchange of integrals ([T] App. A). Mixed $\psi$–$b$ compositions need this checked.
5. **Fixed-rail coupling. Now has a proposal,** Sec. 8: the natural N=2 generalization of [GR] gives a line with uniform $k\to gk$.

## 5. The draft against its sources (Checked)

Script: `scripts/check_literature_inputs.py --dps 25`, 419/419 passing, recorded in `results/data/literature_checks/`. Source formulas are transcribed from the PDFs, each with its equation number. This settles A1.

### 5.1 Component input

| [N2] | Source | Result |
|---|---|---|
| (2) | [PSV](2.9), (2.10); Schwinger–Dyson solved by hand: $2\pi\cot(\pi\delta/2)Jb_\psi^q=1$ | agrees |
| (3) rails, rung | [PSV](2.22)–(2.24), (4.6), (4.7) | agrees, by hand |
| (9), entry by entry, at generic $h$ | [PSV](2.27)–(2.29), (2.33)–(2.35) | agrees, $q=3,5,7$ |
| (15), (16) and charged kernels | [PSV](2.20), (2.21), (2.30), (2.31), (2.36), (2.37), [FGMS](6.2) | agrees, $q=3,5,7$ |
| $H_S$, $H_A$ | [PSV](2.38)–(2.43), all printed roots | agrees |
| (24), $\mathcal C_2$ | [B18](4.9), (4.11), (4.14) | agrees |
| $F_0=z^\delta$; fermion covariance is $-F$ | [PSV](4.9) with $\chi=z/(z-1)$ gives $-z^\delta$; [B18](4.86) | agrees, by hand |

$k_c$ in `check_robin_identities.py`, written from memory, also agrees with [MS](3.73).

### 5.2 The draft's claims of misprints in [B18]

**Both are confirmed.** They are algebraic inconsistencies that we can demonstrate from the paper itself, not an official erratum.

1. **Sign of (4.48), (4.49).** Literal $k^A(1)=-1$, while the text states $k=1$ at $h=1$. And her own (4.51) with (B.8) gives the opposite sign to (4.48): with $(2\Delta-2)\Gamma(2\Delta-2)=\Gamma(2\Delta-1)$ and $1/\Gamma(2\Delta-1)=-\Gamma(2-2\Delta)\sin(2\pi\Delta)/\pi$,
$$\frac{2\Delta+h-1}{2\Delta-2}k^A_{N=0}=+\frac1{\pi^2}\Gamma(-2\Delta)\Gamma(2-2\Delta)\Gamma(2\Delta-h)\Gamma(2\Delta+h)\sin2\pi\Delta\,(\sin2\pi\Delta-\sin\pi h).$$
2. **Sign in (4.39).** From (4.38), (B.4), (B.5) and $\sec\pi(h+1)=-\sec\pi h$:
$$\xi_h=h[(1+\sec\pi h)(F_h-F_{h+1})+(1-\sec\pi h)(F_{1-h}-F_{-h})]=h[(1+\sec\pi h)\varphi_h-(1-\sec\pi h)\varphi_{-h}].$$
(B.4), (B.5) themselves agree with direct integration of (4.26).

**New, not in the draft:** a factor of 2 in (4.90)–(4.92), (4.106), and three problems in App. A (the ordering in (A.9), $Q/2$ in the generators, and $Q=0$). See `literature/notes/bulycheva_2018.md`.

**Normalization of $b$.** [B18](3.9) has $Jb^q=\tan(\pi\Delta)/(4\pi)$, [PSV](2.10) has $1/(2\pi)$. This is not an error. [PSV] has $b_b=\delta b_\psi$, i.e. $G_b=-\partial G_\psi$, whereas [B18]'s $\langle12\rangle$ implies $G_b=-2\partial G_\psi$, i.e. $b_b=2\delta b_\psi$. The Schwinger–Dyson equation fixes $2\pi\cot(\pi\delta/2)(b_b/\delta b_\psi)Jb_\psi^q=1$, which reproduces each.

### 5.3 From [B18](4.106) to [N2](26) (Derived here)

The overlap (4.88) comes from App. G with the literal kernel, $k^A_{\rm literal}=-k_A$. So with the norm (4.72), $\langle\xi,\xi\rangle=-4\pi h\tan\pi h$ at $h=is$, the seed is
$$\frac{\langle F_0,\xi_h\rangle}{\langle\xi_h,\xi_h\rangle}=\frac{\alpha}{2}\frac{k_A(h)}{4\pi h\tan\pi h}.$$
Insert $\xi_h$ with the corrected sign. Send $h\to-h$ in the $\varphi_{-h}$ term, which is allowed on the upward line, using $k_A(-h)=k_S(h)$, and use
$$\frac{1+\sec\pi h}{\tan\pi h}=\cot\frac{\pi h}{2},\qquad\frac{1-\sec\pi h}{\tan\pi h}=-\tan\frac{\pi h}{2}.$$
Resumming each with $1/(1-gk)$ gives
$$\frac{\alpha}{8\pi}\left[\cot\frac{\pi h}{2}\frac{k_A}{1-gk_A}-\tan\frac{\pi h}{2}\frac{k_S}{1-gk_S}\right]\varphi_h.$$
With $\alpha=4\pi\delta/[\nu\tan(\pi\delta/2)]$ (4.89), $\alpha/(8\pi)=\mathcal N_\delta$, **which is [N2](26)**, and the discrete sum gives $D_\delta=\alpha/(4\pi^2)$. The relative sign in (26) depends on the corrected (4.39), and the magnitude needs the $\tfrac12$ of (4.88). Part 4 of the script confirms the seed's leading coefficient is 1 with the $\tfrac12$ and 2 without.

### 5.4 The reduced four-point function at finite $g$ (Checked; key test)

This compares the draft's *auxiliary bulk* side with the *literature's* four-point function, which were computed in entirely different ways.

**Literature side.** [PSV](4.26) with $g$ inserted by the geometric series (our extension; they have $g=1$):
$$f^s=\frac{2\mu(h)k_0^{\psi,s}(h)}{(1-gk_A(h))(1-gk_S(h-1))}\Phi^s_h,\qquad2\mu(h)=\frac{2h-1}{\pi\tan(\pi h/2)}.$$
Closing the contour, the block coefficients are residues. By [T](C.8), a block reduces with $\mathcal I(h)=1/[2\mu(h)\mathcal N(h)]$, while [PSV]'s blocks are $2\mathcal N(h)B_h$, **so $\tan(\pi h/2)$ cancels**:
$$(\mathsf PF_{p=1})(\tau)=-2\sum_{h_*}\mathop{\rm Res}_{h=h_*}\frac{k_0^{\psi,s}(h)}{(1-gk_A(h))(1-gk_S(h-1))}e^{-(h_*-\frac12)|\tau|}.$$
The even sector is the same with $k_0^{\psi,a}$, $k_S(h)$, $k_A(h-1)$, and overall $+2$, using Sec. 7.1.

**Bulk side.** [N2](12), with endpoint $=\tfrac12$ mode: $2b_\psi^2(\mathsf PF_p)=\sigma\ell_\psi w_{11}$, with $w=-\mathrm{Res}\,\mathsf G_{g,p}$.

**Result.** They agree at every root with $\omega<4$ (tolerance $10^{-7}$; typical agreement far better), for $q=3,5$, at $g=0.5,0.9,0.99,1$ for $p=1$ and $g=0.9,0.99,1$ for $p=0$. This includes $\omega=\tfrac12,\tfrac32$ at $g=1$.

**Also**, the block coefficients themselves reproduce [N2](28), including the partner structure $\varphi_H=F_H-F_{H+1}$, at $g=0.7,1$ (part 9):
$$\text{coefficient of }F_H:\ -c_\alpha(H),\qquad\text{coefficient of }F_{H+1}:\ +c_\alpha(H).$$

**What this establishes.** [N2](12), (26), (28) describe the N=2 four-point function of the literature, so the auxiliary system's boundary propagator is the centre-integrated ladder. **Not covered:** entries 12, 21, 22; the $\omega\ge4$ tail; the summed reduction. The even-sector reduction also rests on the prescription of Sec. 7.1.

## 6. The relative cylinder from the super-Casimir (Derived here; checked)

Answers B1, and the speculation that used to be in Sec. 4. Script: `scripts/check_supercasimir.py`, 19/19 passing.

**Assumptions.** $su(1,1|1)$ generators [B18](A.1)–(A.8), corrected as below; a chiral–antichiral pair; $\Delta$ the superfield dimension, $2\Delta=\delta$.

**Step 1: the Casimir, with translations acting first.** As printed, (A.9) has $-L_1L_{-1}$, where $L_1=-\partial_\tau$. That isn't central. The central operator is
$$C=L_0^2-\tfrac14J_0^2-L_{-1}L_1+\tfrac12(G_{+1/2}\bar G_{-1/2}+\bar G_{+1/2}G_{-1/2}).$$
Also, (A.3), (A.5), (A.7) need $Q/2\to Q$ to close with (A.4), and (A.15) can only hold for the pair: $\langle12\rangle^{-2\Delta}$ is invariant under all generators for $Q_1=-2\Delta$, $Q_2=+2\Delta$, the chiral-primary condition of [FGMS]. With these, the implementation reproduces (A.16), (A.17) on all odd components with a residual of $10^{-142}$, so we are computing the right thing.

**Step 2: zero total supermomentum.** Impose that the total $L_1$, $G_{-1/2}$, $\bar G_{-1/2}$ annihilate $f$. Chirality leaves pair coordinates $(\tau_1,\bar\theta_1;\tau_2,\theta_2)$, which is $2|2$. The conditions remove $1|2$. That leaves $1|0$:
$$f=f(x),\qquad x=\langle12\rangle.$$
**So there is no odd relative coordinate.** The relative space is bosonic, and the cylinder is an ordinary $\mathbb R\times S^0$. This is one level below [B18](4.7).

**Step 3: the reduced Casimir.** All terms with a translation on the right drop, and $J_0f=0$. So
$$C\,f(x)=L_0^2f(x),\qquad L_0f(x)=-(x\partial_x+2\Delta)f(x).$$
**Step 4: conjugate.** With $x=e^\tau n$, $n=\pm1$, and $Wf=e^{2\Delta\tau}f$:
$$WCW^{-1}=\partial_\tau^2.$$
On $e^{-H\tau}Y_p(n)$ the eigenvalue is $H^2$, so $\omega=H$, with shadow $-H$.

**Comparison with Thong at $d=1$:**

| | Thong | N=2 |
|---|---|---|
| Reduced Casimir | $L_0^2+L_0$ | $L_0^2$ |
| Weight $a$ | $2\Delta_\phi-\tfrac12$ | $2\Delta=\delta$ |
| $\omega$ | $h-\tfrac12$ | $H$ |
| $WCW^{-1}$ | $\partial_\tau^2-\tfrac14$ | $\partial_\tau^2$ |

**Pattern (literature eigenvalues; our reading).** $C\to L_0^2+(1-\mathcal N/2)L_0$, so $\omega=h-d_s/2$ with $d_s=1-\mathcal N/2$: $h(h-1)$ for $\mathcal N=0$, $h(h-\tfrac12)$ for $\mathcal N=1$ ([MSW](2.24)), $h^2$ for $\mathcal N=2$. This confirms the earlier speculation $d_{\rm eff}=0$. $\mathcal N=2$ is special, as chirality removes the odd relative direction; $\mathcal N=1$ has one (Sec. 9).

**Components.** Expanding, $f(\langle12\rangle)=f(\tau_{12})-2\bar\theta_1\theta_2f'(\tau_{12})+\ldots$. So the zero-supermomentum boson entry is a derivative of the fermion entry. [PSV]'s kernel knows nothing of this, yet its right eigenvector $(1,c)$ on the branch $k_A$ (sector a) or $k_S$ (sector s) has exactly
$$c=-(h-\delta)$$
for all tested $h$, $q=3,5,7$ (part 6). So $G_b=-\partial G_\psi$ extends to all superprimary eigenfunctions, in [PSV]'s normalization. The constant is fixed at $h=0$, where the eigenfunction is the propagator.

**Why the draft has a matrix (Derived here; interpretive).** Ordinary centre integration imposes only $P=0$, so it sees all of $(H,H+\tfrac12,H+\tfrac12,H+1)$: the neutral sector has both bottom and top components, hence two branches and a $2\times2$ matrix with frequency-dependent eigenvectors. Super-centre integration would also impose $G_{-1/2}=\bar G_{-1/2}=0$, which keeps one function $f(x)$ and one scalar kernel eigenvalue, $k_\alpha(H)$. The component weights fit this:
$$a_\psi=\delta-\tfrac12,\quad a_b=\delta+\tfrac12,\qquad\omega_{\rm comp}=H\mp\tfrac12,$$
each $\tfrac12$ from the super value. This is exactly the shift in [N2](16), (19).

**Consequence for C3 (conjecture).** The natural N=2 analogue of Thong's result is one kernel on one bosonic cylinder,
$$\frac{k_\alpha(H)}{1-gk_\alpha(H)}=\sigma L_q\,G^F_{g,p}(H),$$
which is [N2](22), (23): **a single Dirac response, with an algebraic boundary condition $b_F=g\sigma L_qa_F$. The two scalars are its components** in the bosonic reduction, tied by (19), (20). Note that $k_A(H)=k_S(-H)$ isn't even in $\omega$, which is natural for a first-order system and impossible for one Thong scalar. What's missing is the super-centre integral of a superblock.

**Checks.** $H\to-H$ exchanges $A$ and $S$ ([B18](4.50)), on the cylinder $\tau\to-\tau$. The reduced norm is right: $d\chi/[\chi(1-\chi)]\to d\chi/\chi$ at small $\chi$, matching $d\tau$ on the cylinder.

**Not established.** That the super-centre integral exists with the right convergence (B2, B3).

## 7. Exceptional modes and the resolution of C2 (Derived here; checked)

### 7.1 Centre integral in both sectors

Odd sector: [T](C.8), $\mathcal I(h)=\frac{\pi\tan(\pi h/2)}{2\omega\mathcal N(h)}$.

Even sector, new: take [PSV]'s $\Phi^a$ with the same prescription, $B^a_h=B_h$ on $0<\chi<1$, $\Phi^a_h/[2\mathcal N]$ on $\chi>1$, and $-B_h(\chi/(\chi-1))$ on $\chi<0$. Then
$$(\mathsf PB^a_h)(\tau)=\mathcal I_a(h)e^{-\omega|\tau|},\qquad\mathcal I_a(h)=\frac{\pi\cot(\pi h/2)}{2\omega\mathcal N(h)}.$$
**Found and tested by direct quadrature** at three $(h,\tau)$, with the odd sector as control (part 8). Not derived analytically. My first guess had the wrong sign, which quadrature fixed. It requires $\Re h>\tfrac12$, as in [T].

| | $h=1$ | $h=2$ |
|---|---|---|
| $\mathcal I$ (odd) | pole | 0 |
| $\mathcal I_a$ (even) | 0 | pole |

### 7.2 The exceptional families are step functions

$\Phi^s_h$ has poles at odd $h$ and $\Phi^a_h$ at even $h$, **and only on $\chi>1$**. From [PSV](3.6), (3.7), with $x=(\chi-2)/\chi$:
$$\mathop{\rm Res}_{h=1}\Phi^s_h=-2P_0(x),\quad\mathop{\rm Res}_{h=3}\Phi^s_h=-2P_2(x),\quad\mathop{\rm Res}_{h=2}\Phi^a_h=+2P_1(x),\quad\mathop{\rm Res}_{h=4}\Phi^a_h=+2P_3(x),$$
checked (part 8). Meanwhile the measure factors $2\mu(h)$ and $2\mu_a(h)=(2h-1)/[\pi\tan(\pi(h-1)/2)]$ have zeros there, as noted for $\Psi_h$ in [MS] below (3.87). So the product is finite:
$$2\mu(h)\Phi^s_h\to\theta(\chi>1)\quad(h\to1),\qquad2\mu_a(h)\Phi^a_h\to-3P_1(x)\theta(\chi>1)\quad(h\to2).$$

**Lemma.** With $\chi=1/(u^2-s^2)$, $s=\sinh(\tau/2)$,
$$\mathsf P[\theta(\chi>1)P_n(x)]=\frac{2}{2n+1}e^{-(n+\frac12)|\tau|}.$$
*Proof.* $\chi>1$ is $s<|u|<c$, $c=\cosh(\tau/2)$, with $x=1-2(u^2-s^2)$ running over $[-1,1]$. Changing variables,
$$\mathsf P[\ldots]=\int_{-1}^1\frac{P_n(x)\,dx}{\sqrt{2\cosh\tau-2x}}.$$
With $t=e^{-|\tau|}$, $2\cosh\tau-2x=(1-2xt+t^2)/t$, and the Legendre generating function gives $(2\cosh\tau-2x)^{-1/2}=\sqrt t\sum_mt^mP_m(x)$. Orthogonality finishes it. $\square$ Checked numerically at $n\le4$.

So at $h=1$ the interval width is the whole story: $2(c-s)=2e^{-|\tau|/2}$. Combining with the limits above reproduces both exceptional cases of Sec. 5.4 by hand: $-\mathrm{Res}[k_0/D]\cdot2e^{-|\tau|/2}$, and $+3\,\mathrm{Res}[k_0/D]\cdot\tfrac23e^{-3|\tau|/2}$.

### 7.3 Resolution of C2

| Mode | $p$ | $h$ | Coefficient of $B_h$ as $g\to1$ | $\mathcal I$ | Reduced weight |
|---|---|---|---|---|---|
| $A$ | 1 | 1 | $\to0$ | $\infty$ | finite, nonzero |
| $A$ | 0 | 2 | $\to0$ | $\infty$ | finite, nonzero |
| $S$ | 1 | 2 | $\to\infty$ | 0 | finite, nonzero |
| $S$ | 0 | 1 | $\to\infty$ | 0 | finite, nonzero |

- **The $S$ rows are Thong's $h=2$ mechanism** [T](C.10)–(C.12): an enhanced mode that the centre integral kills, leaving a finite remainder.
- **The $A$ rows are its mirror, and were the tension.** At $g=1$ the block expansion on $0<\chi<1$ has no $A$ contribution, so the draft's $c_A\propto\cot(\pi H/2)\to0$ is correct. **But the four-point function does contain the mode.** Its contribution is a finite step function, supported only on $\chi>1$:
$$F_{p=1}\supset-\mathop{\rm Res}_{h=1}\left[\frac{k_0^{\psi,s}}{D^a}\right]\theta(\chi>1).$$
The draft's reconstruction lives on $0<z<1$, i.e. $\chi<0$, which can't see it. That is exactly the loophole in its phrase "in this particular seed and ordering".

So **[T](3.13) isn't violated**, and the draft's two statements are compatible. The mode has zero OPE coefficient in the usual sense, yet a finite reduced weight, because the family isn't a function on $0<\chi<1$ at all. It only appears in the $ijij$ ordering of [MS](3.57).

**Physics ([FGMS] Sec. VI, stated there as an expectation).** The $S$ multiplet "should lead to the superschwarzian", and the $A$ multiplet "should correspond to operators in the IR theory", from the extra symmetry (6.6) with $f\ne h$. Our result fits this: no divergence and a finite contribution. It adds its explicit form in the four-point function, which they leave open.

**Caveat.** At $g=1$ exactly, $\mathcal I=\infty$ and the coefficient is 0. So the family must be defined as the limit $g\to1^-$, or directly as the step function. The physical interpretation of a contribution in only one ordering has not been examined.

## 8. A supersymmetric line for $g$ (Derived here at leading large $N$)

Addresses D2. The draft says "no microscopic supersymmetric line parametrized by $g$ is assumed", but [GR]'s construction appears to carry over.

**Proposal.** Replace the kinetic term by a superspace bilocal,
$$S_0=\int d\mu_1d\mu_2\,\bar\Psi(1)K(\langle12\rangle)\Psi(2),$$
with $K$ the inverse of $\bar G=\bar b\,\mathrm{sgn}\langle12\rangle/|\langle12\rangle|^{\delta}$, keeping the superpotential. In components this should give UV propagators $\bar G_\psi$ and $\bar G_b=-\partial\bar G_\psi$. I haven't carried out this reduction explicitly.

**Schwinger–Dyson.** $G_i^{-1}=\bar G_i^{-1}-\Sigma_i$, with $\Sigma$ as in [PSV](2.5), (2.6). Try the conformal form, with $b_\psi/\bar b_\psi=b_b/\bar b_b=r$. In frequency space $G_i\bar G_i^{-1}=r$, and $G_i\Sigma_i=-g$ with
$$g\equiv\frac{Jb_\psi^{q-1}b_b}{(Jb_\psi^{q-1}b_b)_{\rm IR}},$$
**the same for both equations**, since the IR solution makes each product $-1$. So
$$g=1-r.$$
This is [T](C.9) with $r=2\bar b$. $r=1$ is free, and $r\to0$ is the IR.

**Kernel.** All three entries scale the same way: $K^{11}\propto Jb_\psi^{q-1}b_b$, $K^{12}\propto Jb_\psi^q$, $K^{21}\propto Jb_b^2b_\psi^{q-2}$, with $b_b/b_\psi=\delta$ fixed. So
$$\widehat{\mathsf k}\to g\widehat{\mathsf k},$$
which is **exactly the draft's formal $g$** in [N2](10), with fixed rails.

**Consequences (open).**

- **The $H=1$ multiplets lift for $g<1$**, as $k_\alpha(1)=1\neq1/g$. That's consistent with a nonlocal theory with no local stress tensor, and it's the N=2 version of [GR]'s $h_0$.
- **$g_*$ may be physical.** On this line, the lowest $p=0$ operator would reach $h=\tfrac12$ at $g_*<1$. Near $g_*$ the poles at $h$ and $1-h$ would pinch the contour at $\Re h=\tfrac12$, which fits the $1/\omega$ growth of the residue in Sec. 3.3. In [GR] the analogous $h$ stays above $\tfrac12$. So the N=2 line may split into two segments. **This is a place where $\mathcal N=2$ differs from $\mathcal N=0$, and it needs its own study.**
- The seed's pole at $h=\delta<\tfrac12$ says the UV is a generalized free field that is non-normalizable in that channel, as in [B18](4.103).
- **Not checked:** the component reduction of $S_0$, reflection positivity, or anything at order $1/N$.

## 9. Stepping stones (Checked in part)

**Complex SYK at zero charge.** $k^S(h)\mathcal R_0(h-\tfrac12)=4^{-\nu}$ for $\nu=1-2/q$, tested at $q=4,6,8,12$ and complex $h$ (part 10), against $k^A\mathcal R_1=(q-1)4^{-\nu}$ ([T]). So a **single scalar with a parity-dependent Robin coefficient**, $\kappa_{\rm even}=g4^{-\nu}$, $\kappa_{\rm odd}=g(q-1)4^{-\nu}$. It has $h=1$ for the $U(1)$ charge, and $a+J=2/q-\tfrac12<0$ for $q>4$. So it's the simplest system with the $a+J<0$ problem and an $h=1$ mode, with no mixing. Thong's "beyond the single scalar" (p. 6) isn't needed at zero charge, though the parity dependence couples the two endpoints. The seed normalization and block reduction in this channel haven't been carried out; Sec. 7 should apply.

**N=1 has an odd relative coordinate** (Sec. 6, `literature/notes/murugan_stanford_witten_2017.md`). So it isn't a sub-case of the N=2 structure, and tests something different.
