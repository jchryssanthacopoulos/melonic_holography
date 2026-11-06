# Gu, Kitaev, Sachdev, Tarnopolsky 2019

**Citation.** Y. Gu, A. Kitaev, S. Sachdev, G. Tarnopolsky, *Notes on the complex Sachdev-Ye-Kitaev model*, JHEP 02 (2020) 157, arXiv:1910.14099v3. PDF: `papers/gu_kitaev_sachdev_tarnopolsky_2019.pdf`. Tag: **[GKST]**. The journal reference is from memory, and hasn't been checked.

**Read: targeted only.** The abstract, contents, Sec. 5.1 and App. B. Everything else in this 75-page paper was not read. **Do not treat this note as a full review.**

## Main question

Properties of complex SYK with a $U(1)$ charge: the charge's relation to the IR asymmetry, compressibility, density of states, and a bulk picture of the entropy.

## Physical system

Complex SYK at arbitrary charge, with spectral asymmetry $\theta$ or $\mathcal E$.

## Main analytical results used here

- **Spectrum at $\theta=0$** (B.1)–(B.3): a symmetric matrix $W_\Sigma$ with eigenvectors $(1,\pm1)$, giving $k^A$ and $k^S$. $k^A$ is the SYK kernel, with $h_0^A=2$. $k^S$ has $h_0^S=1$, the $U(1)$ charge.
- **At $\theta\neq0$ the channels mix** (B.4), (B.5), via a $2\times2$ kernel with a square root. At integer $h$ the eigenvalues don't depend on $\theta$, so $k^A(2,\theta)=k^S(1,\theta)=1$ for all $\theta$: **the energy and charge dimensions are protected**.
- **There's a critical value**, $\cos2\theta_c=1-2\Delta$ (B.6), above which $h^S_{n\ge1}$ disappear.
- **Bulk picture** (Sec. 5): the zero-temperature entropy equals that of a free massive Dirac (or Majorana) fermion on $H^2$, with $M=\tfrac12-\Delta$, in an imaginary electric field.

## Limitations

For our purposes, those of a partial reading.

## Relationship to our project

- **Check**: (B.3) agrees with [B18](B.9), and $k^S\mathcal R_0=4^{-\nu}$ (part 10 of `scripts/check_literature_inputs.py`).
- **D3 (nonzero charge).** (B.4), (B.5) are a worked example of $\theta$-dependent channel mixing, the structure a charged Thong construction would have to reproduce. This suggests a matrix Robin condition depending on $\theta$. Speculative.
- **Do not conflate the Dirac fields.** $M=\tfrac12-\Delta$ here describes a *single-site* fermion of dimension $\Delta$, with boundary weights $\Delta$ and $1-\Delta$. The draft's Dirac field has $\mu=\tfrac12-\delta$, with weights $\delta$ and $1-\delta$, and represents a *pair* with centre integrated out. Similar equations, different objects.
