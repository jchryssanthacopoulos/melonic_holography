# Ishibashi, Wald 2004

**Citation.** A. Ishibashi, R. M. Wald, *Dynamics in non-globally-hyperbolic static spacetimes III: Anti-de Sitter spacetime*, Class. Quant. Grav. 21 (2004) 2981, arXiv:hep-th/0402184v4. PDF: `papers/ishibashi_wald_2004.pdf`. Tag: **[IW]**.

**Read:** the abstract, the setup in Sec. 1, and Sec. 3.3 through Theorem 3.2. Secs. 2 and 3.1–3.2 were only searched; Sec. 4 was not read.

## Main question

AdS isn't globally hyperbolic, so what are the possible dynamics of a field on it?

## Physical system

Scalar, electromagnetic and gravitational perturbations on $(n+2)$-dimensional AdS, $n\ge1$; the text says the results extend straightforwardly to $n=0$ for scalars. So AdS$_2$ is covered, though it's not the main case.

## Important definitions

Dynamics is set by a positive self-adjoint extension $A_E$ of the spatial operator $A$. Radial problem: index $\nu$, with $\sigma$ set by the harmonic, with hypergeometric solutions (144), (145).

## Main analytical results

- $\nu^2\ge1$: unique extension.
- $0<\nu<1$: a one-parameter family, labelled by $b_\nu/a_\nu$, the ratio of the two boundary falloffs. These are "generalized Robin" conditions, and Dirichlet is the Friedrichs extension.
- $\nu=\tfrac12$ gives ordinary Robin (168). They note the double-trace connection in footnote 5.
- **Theorem 3.2, Eq. (169)**: positive if and only if $b_\nu/a_\nu\ge-\frac{\Gamma(-\nu)}{\Gamma(\nu)}\frac{\Gamma(\zeta^0_{\nu,\sigma})^2}{\Gamma(\zeta^0_{-\nu,\sigma})^2}$.

The form of (169) is Thong's condition $\kappa>\mathcal R_J(0)$, with his $\mathcal R_J(0)=\frac{\Gamma(-\nu)}{\Gamma(\nu)}[\Gamma(\frac{J+\Delta_+}{2})/\Gamma(\frac{J+\Delta_-}{2})]^2$. The overall sign in the printed extraction of (169) needs checking against the PDF, and I haven't mapped their sign convention for $b_\nu$ or variables $\sigma,\zeta$ in detail. **Treat this as matching in form only.**

## Numerical methods

None.

## Limitations

Single-channel: one field, one boundary parameter. No matrix boundary conditions. Assumes a static, real problem.

## Relationship to our project

- Foundation for Thong's stability statement.
- **It doesn't cover the draft's system.** The matrix condition couples two fields. That's why `docs/derivations.md` 3.1, 3.2 prove self-adjointness and $\omega^2\ge0$ directly. Those remain project-derived.
- The sector $\nu=\delta$, $\rho=\delta$ with standard quantization is within $0<\nu<1$, so each scalar separately has a boundary-condition family.
