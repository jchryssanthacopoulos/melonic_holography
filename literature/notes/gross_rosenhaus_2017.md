# Gross, Rosenhaus 2017

**Citation.** D. J. Gross, V. Rosenhaus, *A line of CFTs: from generalized free fields to SYK*, JHEP 07 (2017) 086, arXiv:1706.07015v2. PDF: `papers/gross_rosenhaus_2017.pdf`. Tag: **[GR]**.

**Read:** Secs. 1, 2 (through the free energy discussion), and 3.1 through (3.35). Sec. 4 and App. A were not read.

## Main question

Is there a variant of SYK that is $SL(2,\mathbb R)$ invariant at every coupling, so that weak and strong coupling can be connected?

## Physical system

Conformal SYK (cSYK): replace the kinetic term by a bilocal one that already gives $\chi$ dimension $\Delta=1/q$ (1.2), (3.7), keeping the SYK interaction (3.9). $J$ is then dimensionless, and there's a line of fixed points at least at leading $1/N$.

## Main assumptions

Large $N$. UV finiteness is argued from fermion antisymmetry and randomness; conformal invariance to all orders in $1/N$ is left as possible.

## Main analytical results

- **Propagator**: $G=b\,\mathrm{sgn}\tau/|\tau|^{2\Delta}$ with $b^q/(1-2b)=\frac{1}{2\pi J^2}(1-2\Delta)\tan\pi\Delta$ (3.12). $b\to\tfrac12$ at weak coupling, with $b\to0$ at strong.
- **Kernel** (3.25): $\tilde g(h)=(1-2b)g(h)$, where $g$ is the SYK eigenvalue. So $k\to(1-2b)k$ uniformly.
- **OPE** (3.31): the SYK answer, with $(1-2b)^{-2}$.
- **Spectrum**: $h\to2\Delta+2n+1$ as $J\to0$ (3.27). The lowest $h_0\to2$ as $J\to\infty$ with $c_0^2\to\infty$ (3.28)–(3.33). The text says $h_0$ "increases from 1", which isn't the $n=0$ value of (3.27); we go with (3.27), as Thong does.
- **Bulk**: a bilocal action is equivalent to a free scalar with a boundary field (2.5), or a continuous tower (2.9). The dual of cSYK is non-gravitational, on rigid AdS$_2$ (3.34).
- The $O_0$ correlators are conformal in cSYK, and break it in SYK.

## Numerical methods

Solving (3.12) for $b(J)$, figure 1.

## Limitations

Leading $1/N$; the bilocal action needs a UV origin (auxiliary tower); no time-reparameterization symmetry, by construction, and no local stress tensor.

## Relationship to our project

- **Meaning of $g$**: this is what makes Thong's $g$ microscopic, with $g=1-2\bar b$ [T](C.9).
- **D2 (is $g\neq1$ meaningful in N=2?).** The draft assumes no such line. The argument here extends directly to N=2, since all three kernel entries scale as $Jb^q$. **So a superspace bilocal action gives a line with $g=1-b/\bar b$** at the level of the Schwinger–Dyson equations. Derived in `docs/derivations.md` Sec. 8. This makes the draft's formal $g$, and also $g_*$, potentially physical.
- **Comparison**: [GR] has a tower of bulk fields (3.34), Thong has one field. They describe the same four-point data for the same theory.
