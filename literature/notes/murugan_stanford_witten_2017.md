# Murugan, Stanford, Witten 2017

**Citation.** J. Murugan, D. Stanford, E. Witten, *More on supersymmetric and 2d analogs of the SYK model*, JHEP 08 (2017) 146, arXiv:1706.05362v1. PDF: `papers/murugan_stanford_witten_2017.pdf`. Tag: **[MSW]**. The journal reference is from memory, and hasn't been checked.

**Read: targeted only.** The abstract, contents, Sec. 2 (the Casimir and its eigenvalue), and the start of Sec. 3.6 (the N=1 kernel). The paper is 106 pages, and Secs. 4–9 on the 2d models were not read. **Do not treat this note as a full review.**

## Main question

Give a basis of (super)conformal eigenfunctions for four-point functions, clarify 1d N=1 SYK, and build 2d SYK analogues.

## Physical system

1d N=1 with $\hat q$-fold superspace interactions, and 2d bosonic and supersymmetric models.

## Important definitions

- $osp(1|2)$ generators on $(t,\theta)$ (2.2), and Casimir (2.3), with one-particle value $\Delta^2-\Delta/2$ (2.4).
- Two-particle Casimir (2.5). **Eigenvalue $h(h-\tfrac12)$**, invariant under $h\to\tfrac12-h$ (2.24).
- N=1 has a pair of functions $F,G$ of $\chi$ with $G=\chi^hF$ or similar, each with an ordinary Casimir (2.25)–(2.27).
- Note: their $L_{-1}=-\partial_t$ is translation, where [B18] calls it $L_1$.

## Main analytical results used here

- N=1 kernel (3.49)–(3.52), and eigenfunctions of bosonic and fermionic type (3.53).
- $k^F(\tfrac32)=k^B(-1)=1$ (near line 1154 of the extraction), the supercurrent multiplet.

## Limitations

For our purposes, those of a partial reading.

## Relationship to our project

**D1: the pattern across $\mathcal N$ (derived here).** Reordering the Casimir with translations on the right, and dropping them at zero total supermomentum, gives
$$C\to L_0^2+\left(1-\frac{\mathcal N}{2}\right)L_0,$$
with $L_0\to-h$. So:

| $\mathcal N$ | Eigenvalue | Shadow | Source | $\omega$ | $d_s$ |
|---|---|---|---|---|---|
| 0 | $h(h-1)$ | $1-h$ | [MS](3.60) | $h-\tfrac12$ | 1 |
| 1 | $h(h-\tfrac12)$ | $\tfrac12-h$ | [MSW](2.24) | $h-\tfrac14$ | $\tfrac12$ |
| 2 | $h^2$ | $-h$ | [B18](4.10) | $h$ | 0 |

Here $d_s=1-\mathcal N/2$ plays the role of $d$ in $\omega=h-d/2$. The eigenvalues are all from the literature; **$\omega$ and $d_s$ are the project's reading**. The full N=2 derivation is in `docs/derivations.md` Sec. 6. The N=1 case above was done by hand from (2.3) and has not been run through code.

**N=1 is structurally different.** Its pair coordinates are $2|2$ and supertranslations remove $1|1$, leaving $1|1$. So **N=1 does have an odd relative coordinate**, a true super-cylinder, whereas chirality removes it in N=2. N=1 isn't simply a simpler N=2.
