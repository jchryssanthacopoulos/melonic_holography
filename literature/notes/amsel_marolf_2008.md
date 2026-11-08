# Amsel, Marolf 2008

**Citation.** A. J. Amsel, D. Marolf, *Supersymmetric Multi-trace Boundary Conditions in AdS*, Class. Quant. Grav. 26 (2009) 025010, arXiv:0808.2184v2. Tag: **[AM]**. The PDF isn't in `papers/`; it was read from an arXiv download, and it would be worth adding.

**Read:** Sec. I, the $d=4$ analysis in Sec. III, and Sec. VI. Secs. II, IV, V and the appendices were not read closely.

## Main question

Which multi-trace boundary conditions for bulk scalar multiplets preserve supersymmetry in AdS$_d$?

## Physical system

Free massive Dirac fermions and scalars on fixed AdS$_d$, $d\ge2$ for fermions, and $d=3,4,5$ for the supersymmetry analysis. Table I gives the field content.

## Main assumptions

Standard inner products. Free fields, aside from possibly non-linear boundary conditions.

## Main analytical results

- Fermions with $|m|<\tfrac12$ have a choice of boundary conditions. For $|m|\ge\tfrac12$ the condition is unique.
- Scalar and fermion masses in a multiplet: $m_{\phi,\pm}^2=m_\psi^2\pm m_\psi-d(d-2)/4$ (6.1), said to hold in $d=2$ too.
- Boundary superfields (3.24), in $d=4$:
$$\Phi_-=\alpha_A+\theta a_\psi+\tfrac12\theta\theta(1-2m)\beta_B,\qquad\Phi_+=\alpha_B-\theta\beta_\psi-\tfrac12\theta\theta(1+2m)\beta_A.$$
- Boundary conditions $\Phi_-=\delta W(\Phi_+)/\delta\Phi_+$ (3.27) conserve the inner product and preserve boundary Poincaré supersymmetry, half the supercharges. They correspond to adding $\int d^3S\,d^2\theta\,W$ to a dual CFT.
- The linear example is $\Phi_-=q\Phi_+$ (3.32).
- **"Supersymmetric boundary conditions always relate two distinct bulk scalar fields"**, so single-scalar multi-trace conditions have no supersymmetric generalization (Sec. I).
- **"In general, linear boundary conditions can preserve only the Poincaré subalgebra of supercharges"** (Sec. VI). Full superconformal symmetry needs special masses and special $W$.

## Numerical methods

None.

## Limitations

No $d=2$ supersymmetry analysis (they cite earlier work for single-trace AdS$_2$). $\mathcal N=1$ only. Boundary Poincaré symmetry is the preserved subalgebra they consider.

## Relationship to our project

This is the closest precedent for Secs. 5, 6 of `research/pdf/n2_dirac_ladders.pdf`, and I'd only found it in the literature search for the paper.

- **The structure is the same.** Our $\Phi_\mp$ are the $d=2$ analogue of their superfields, which they construct for $d=3,4,5$, with $(1\mp2m)\to2\delta,2\nu$ at $m=\mu$. Our masses are their (6.1) at $d=2$.
- **The difference.** Their condition has constant $q$ in flat space, which breaks dilatations and preserves the Poincaré half. Ours has constant $c$ on the cylinder, which preserves dilatations and breaks translations, and it keeps a different half: $\{G_{-1/2},\bar G_{+1/2}\}$, closing on $2L_0-J_0$ rather than $P$.
- **Consequence for the draft.** The matrix condition [N2](11), which looked specific to this problem, is a known structure in a different frame. The Weyl shift in the top components accounts for the entry $2\mu$ of $\mathsf M_p$.
- **Consistent with Thong's remark** that a single scalar can't be the supersymmetric statement.
