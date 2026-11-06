# Bulycheva 2017

**Citation.** K. Bulycheva, *A note on the SYK model with complex fermions*, JHEP 12 (2017) 069, arXiv:1706.07411v1. PDF: `papers/bulycheva_2017.pdf`. Tag: **[B17]**.

**Read:** Sec. 1, parts of Sec. 4, and Sec. 7. Secs. 2, 3, 5 and 6 were searched only. These are mostly superseded for our purpose by [B18].

## Main question

Do the conformal four-point analysis for complex fermions.

## Physical system

Complex SYK at zero charge. Same $\Delta=1/q$. Two channels, by $\chi\to\chi/(\chi-1)$ parity.

## Main analytical results

- Symmetric-channel eigenvalue $k^S$, with $k^S(h,\tfrac14)=\tfrac12\cot(\tfrac\pi2(h-\tfrac12))/(h-\tfrac12)$ at $q=4$ (4.29).
- **$U(1)$ mode**: $k^S=1$ at $h=0$ (equivalently $h=1$), for all $q>2$. Operators $O_k^S=\bar\psi\partial^{2k}\psi$ with $h\approx2k+2\Delta$ (4.31)–(4.33).
- $T$-odd functions exist but don't contribute at leading order.
- The $U(1)$ mode is an eigenfunction of the retarded kernel, with zero Lyapunov exponent.
- Sec. 7: the charge mode gives a divergence that, unlike $h=2$, is not cured by weakly breaking the symmetry, since $U(1)$ is exact.

## Limitations

Zero charge only. Later work ([GKST]) revisits the charge mode with care.

## Relationship to our project

**A new result for D1 (derived here, checked).** Complex SYK's symmetric channel is exactly Thong's $J=0$ sector of the same scalar:
$$k^S(h)\mathcal R_0(h-\tfrac12)=4^{-\nu},\qquad\text{compare }k_c(h)\mathcal R_1(h-\tfrac12)=(q-1)4^{-\nu}\ [T].$$
Tested for $q=4,6,8,12$ in part 10 of `scripts/check_literature_inputs.py`. Consequences:

1. **Single scalar, parity-dependent Robin coefficient**, at zero charge. So Thong's remark that complex SYK goes "beyond the single scalar" isn't needed at this level. The even and odd parities need different $\kappa$, which on the strip is a coupling between the two endpoints.
2. **It has $a+J=2/q-\tfrac12<0$ for $q>4$**, the same failure as the N=2 $p=0$ sector, and a mode at $h=1$. So it's the cleanest test case for continuing past Thong's validity range, with no mixing involved.
