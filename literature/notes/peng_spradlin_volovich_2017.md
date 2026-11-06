# Peng, Spradlin, Volovich 2017

**Citation.** C. Peng, M. Spradlin, A. Volovich, *Correlators in the N = 2 Supersymmetric SYK Model*, JHEP 10 (2017) 202, arXiv:1706.06078v2. PDF: `papers/peng_spradlin_volovich_2017.pdf`. Tag: **[PSV]**.

**Read:** in full.

## Main question

Compute the leading $1/N$ four-point functions of N=2 SYK in components, and check they are consistent with supersymmetry.

## Physical system

$N$ complex fermions $\psi_i$ and auxiliary complex bosons $b_i$, odd $q$ (2.1), with $\langle C\bar C\rangle=(q-1)!J/N^{q-1}$ (2.2). Note the coupling is $J$, not $J^2$.

## Definitions and solution

(2.9), (2.10): $G^\psi=b_\psi\,\mathrm{sgn}\,\tau/|\tau|^{2\Delta_\psi}$, $G^b=b_b/|\tau|^{2\Delta_b}$, with
$$\Delta_\psi=\frac1{2q},\quad\Delta_b=\frac{1+q}{2q},\quad b_\psi^q=\frac{\tan(\pi/2q)}{2\pi J},\quad b_b=\frac{b_\psi}{q}.$$
(2.14): $\alpha=\frac{q-1}{q}\frac{\tan(\pi/2q)}{2\pi}$.

**Label trap.** A 2-component eigenfunction is called "symmetric" or "antisymmetric" by its *fermion* entry. But the four-point functions in Sec. 4 are labelled by the symmetry of $\Phi_h$ under $\chi\to\chi/(\chi-1)$, which is opposite. So $f^s$ (with $\Phi^s$) has denominators in $k^{a,\pm}$: Eq. (4.26). In the project's notation $p=1$ is PSV's "a" for kernels and "s" for $\Phi$. See `docs/derivations.md` Sec. 1.

## Main results

- **Diagonal kernels** (2.12)–(2.21), for $\psi b$ channels, labelled $k^{s,d},k^{a,d}$.
- **Matrix kernels** (2.22)–(2.25), with entries (2.27)–(2.29) and (2.33)–(2.35). $K^{22}=0$.
- **Eigenvalues** (2.30), (2.36), related by (2.31), (2.37): $k^{a,+}(h-\tfrac12)=k^{a,d}(h)$, $k^{a,-}(h+\tfrac12)=k^{s,d}(h)$, and $s\leftrightarrow a$. Footnote 1: $\tilde k^\pm$ from the square-root formula are the same *set* as $k^\pm$, but not pointwise.
- **Spectra at $q=3$** (2.38)–(2.43), organized into $\{h-\tfrac12,h,h+\tfrac12\}$ towers. Lowest: $(1,\tfrac32,2)$ in both.
- **Chaos** (2.54), (2.55): $k_R^-(-1)=1$, maximal.
- **Eigenfunctions** (3.6)–(3.11), with $\Phi^a$ new relative to Maldacena–Stanford. Inner product (3.12), the $N=0$ measure $d\chi/\chi^2$. Discrete series: even $h$ for $\Phi^s$, odd $h$ for $\Phi^a$ (3.13), (3.14).
- **Four-point functions.** Seeds (4.6)–(4.9). Overlaps (4.10)–(4.13). Resummation (4.25), (4.26), results (4.27)–(4.31), (4.35)–(4.39). The $\psi b\bar\psi\bar b$ function is obtained as $\partial_{\tau_4}$ of an auxiliary "four-fermion" function, (4.42), (4.56).
- **Divergences** (end of Sec. 4.4): come only from $(h^{s,+}_0,h^{s,d}_0,h^{a,-}_0)$. So only one of the two $(1,\tfrac32,2)$ multiplets diverges.

## Numerical methods

Kernels are analytic. Spectra are root-finding (figure 3).

## Limitations

The enhanced multiplet is left out ("we will not do so here"). Sec. 5 on six points is only a sketch.

## Relationship to our project

This is the draft's source for its component input. **The comparison is now done, and the draft agrees** (`scripts/check_literature_inputs.py`, parts 1, 2, 7, 9; `docs/derivations.md` Sec. 5):

| Draft | PSV | Result |
|---|---|---|
| [N2](2) | (2.9), (2.10) | agrees (by hand) |
| [N2](3) kernels | (2.22)–(2.24), (4.6)–(4.7) | agrees (by hand) |
| [N2](9), entrywise | (2.27)–(2.29), (2.33)–(2.35) | agrees, $q=3,5,7$ |
| [N2](15), (16) | (2.20), (2.21), (2.30), (2.36) | agrees |
| $q=3$ roots | (2.38)–(2.43) | agrees |
| [N2](28) at finite $g$ | (4.10), (4.11), (4.26) | agrees, $q=3,5$ |
| "$-F$" with $\chi=z/(z-1)$ | (4.9) gives $F_0=-z^\delta$ | agrees (by hand) |

PSV's formulas with $g$ inserted are our extension: they treat only $g=1$.
