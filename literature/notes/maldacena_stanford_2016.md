# Maldacena, Stanford 2016

**Citation.** J. Maldacena, D. Stanford, *Remarks on the Sachdev-Ye-Kitaev model*, Phys. Rev. D 94 (2016) 106002, arXiv:1604.07818. PDF: `papers/maldacena_standford_2016.pdf` (the filename has a typo). Tag: **[MS]**.

**The PDF is arXiv v1, titled "Comments on the Sachdev-Ye-Kitaev model", and its equation numbers differ from the published version.** Thong cites the published numbering. For example, Thong's [MS](3.35) for the kernel is (3.73) in v1, and Thong's (3.45) for the ladder sum is (3.83)–(3.84) in v1. **All numbers below are v1.** It would be worth pulling the current arXiv version if we need to check further cross-references; the formula used here is checked by its special cases.

**Read:** Secs. 3.2.1–3.2.5. Secs. 1, 2, 3.3 and 4–7 were not read.

## Main question

Analyze the SYK model at large $N$: four-point functions, the soft mode, and the possible bulk.

## Physical system

$N$ Majorana fermions with random $q$-body couplings, $q$ even.

## Important definitions

- $\chi=\tau_{12}\tau_{34}/(\tau_{13}\tau_{24})$, (3.56).
- $\mathcal C=\chi^2(1-\chi)\partial_\chi^2-\chi^2\partial_\chi$, eigenvalue $h(h-1)$ (3.60).
- $\Psi_h$: (3.64), with $\Psi_h=\Psi_{1-h}$; on $0<\chi<1$, (3.65), (3.66). Integral form (3.67). Spectrum $h=\tfrac12+is$ and $h=2n$.
- Seed $\mathcal F_0=-G_{13}G_{24}+G_{14}G_{23}$ (3.41). Kernel (3.44).

## Main analytical results

- **Kernel eigenvalue** (3.73):
$$k_c(h)=-(q-1)\frac{\Gamma(\tfrac32-\tfrac1q)\Gamma(1-\tfrac1q)}{\Gamma(\tfrac12+\tfrac1q)\Gamma(\tfrac1q)}\frac{\Gamma(\tfrac1q+\tfrac h2)\Gamma(\tfrac12+\tfrac1q-\tfrac h2)}{\Gamma(\tfrac32-\tfrac1q-\tfrac h2)\Gamma(1-\tfrac1q+\tfrac h2)},$$
with $k_c(h)=k_c(1-h)$ and special cases (3.75)–(3.77).
- $k_c(0)=-(q-1)$ from the Schwinger–Dyson equation, $k_c(2)=1$ from reparameterizations.
- Overlap $\langle\mathcal F_0,\Psi_h\rangle=\alpha_0k_c(h)/2$ (3.82), and the ladder sum (3.83).
- Near $\chi=1$, $\Psi_h$ has an expansion with $\tan(\pi h/2)$ (3.68). **The poles of $\Psi_h$ at $h=1+2n$ are cancelled by zeros of $1/\tan(\pi h/2)$** (below (3.87)). This is the mechanism behind the $A$-branch result in `docs/derivations.md` Sec. 7.
- Footnote 5: the replacement $1-k_c\to1-gk_c$, with $g<1$, removes the $h=2$ divergence. Thong's $g$ uses this.

## Limitations

Conformal limit only for this section. The $h=2$ sector needs Sec. 3.3, which we haven't read.

## Relationship to our project

- **`kc_syk`** in `scripts/check_robin_identities.py` was written from memory. **It agrees with (3.73)** (part 5 of `scripts/check_literature_inputs.py`), including $q=4$ and the two special values. That to-do item is closed.
- **Template**: PSV and Bulycheva both follow Secs. 3.2.1–3.2.5.
- **Next**: read Sec. 3.3 before working on the physical $H=1$ coefficient (C5).
