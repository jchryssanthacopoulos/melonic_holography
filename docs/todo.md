# To do

Last updated 2026-11-06. Research questions (A1, B1, ...) are in [`research_questions.md`](research_questions.md). Derivation sections are in [`derivations.md`](derivations.md).

## Now: the main remaining gap

- [ ] **B2, C3.** Compute the super-centre integral of a superblock: integrate $\varphi_H(\chi)$ with $\chi$ the super cross ratio over $(Y,\theta_Y,\bar\theta_Y)$, using $\langle ij\rangle$ of [B18](2.12). Expect one exponential $e^{-H|\tau|}$ per parity. Test whether $\mathcal T_2$ of [N2](25) is its inverse. This decides whether one Dirac field is the right statement (derivations 6).
- [ ] **Reduced super-kernel.** Derive $K^{(0)}$ on functions of $\langle12\rangle$ directly from [B18](4.42), the analogue of [T](3.3)–(3.6). Do this first if it's easier: it should give $k_\alpha(H)$ without any component mixing.
- [ ] **Other entries of [N2](12).** Extend `check_finite_g` to 12, 21, 22, using $m^s$, $p^s$, $b^s$ of [PSV](4.26), (4.35). The boson rail needs the centre integral with $\Delta_b$.
- [ ] **Derive $\mathcal I_a(h)$** analytically, following [T](C.6)–(C.8). It's only been found by quadrature.
- [ ] **D5.** Run a literature search for related work, since derivations 6–8 may not be new. Include: N=2 SYK and AdS$_2$ Dirac; super-Casimir at zero momentum; supersymmetric cSYK; "Robin" or "double-trace" in supersymmetric SYK.

## Next: physics of the new results

- [ ] **C4.** Test the draft's boundary domain against the AdS$_2$ Killing spinors ([CGS] (A.7)), not only $E_+$. Find how many supercharges survive. Try extending [CGS] Type I and Type II to $m_F=\mu$.
- [ ] **D2.** Write out the component form of the superspace bilocal $S_0$ (derivations 8). Check supersymmetry and the sign of $\bar G_b$.
- [ ] **D2.** Study $g_*$ on the proposed line: the pinch at $h=\tfrac12$, whether $h$ goes complex or the pole moves through, and the two-segment picture. Check at $g\to0$ against generalized free fields.
- [ ] **C2.** Understand what a contribution in only one ordering means, and how it relates to [FGMS] (6.6). Compute the same step-function term at $h=3,4,\ldots$ for $g\neq1$, where it's generically absent.
- [ ] **C5.** Read [MS] Sec. 3.3 and [FGMS] Sec. IV, V. Derive the physical coefficient of the $S$ multiplet from the super-Schwarzian, and compare with [N2](33).
- [ ] **D1.** Work through complex SYK as a test case: normalization, block reduction, positivity and $g_*$ with no mixing.

## Then: remaining verification

- [ ] **A2.** Test [N2](27) at order $g^2$ against a direct two-rung integral.
- [ ] Test (27) at finite $g$: contour against root sum. The draft reports this; we have not reproduced it.
- [ ] Re-derive [N2](25): derive $\mathcal A$ as an inversion of [T](C.6), since it is not in arXiv v1.
- [ ] Re-derive the rail relation $r_b/c_b=-4\delta^2\ell_\psi$. $r_b$ isn't printed in the draft, so it hasn't been checked.
- [ ] **C1.** Prove or falsify positivity of the residues in general, and $g_*<1$ for all odd $q\ge3$. Currently tested at $q=3,5,7$ only.
- [ ] **B3.** Decide if the mixed continuation is needed in the super formulation.
- [ ] Check [IW](169) against the PDF and map its variables to [T]'s. The note only claims a match in form.

## Literature

- [ ] Finish reading, in this order: [MS] Sec. 3.3; [FGMS] Secs. I–IV; [GR] Sec. 4; [CGS] Secs. 3.3, 4; [MSW] Secs. 2–3 properly; [GKST] Sec. 2–3; [IW] Secs. 3.1, 3.2.
- [ ] Get arXiv's current version of Maldacena–Stanford. The PDF in `papers/` is v1 ("Comments on..."), and Thong cites later numbering.
- [ ] Find the draft's calculation records, `n2_susy_robin_extension.md` and `n2_completion_supplement.md`, and Thong's earlier "supplied draft", if available. Add them to `research/`.

## Project setup

- [ ] Decide on `git init`. The repository isn't under version control, and CLAUDE.md asks for commits at milestones and for code versions in experiment records.
- [ ] Set up a project environment. `.venv` doesn't exist, and SciPy in the current Anaconda install fails to import (NumPy binary mismatch). Current checks need `mpmath` and `sympy`, both working in system `python3`.
- [ ] Add `scripts/memwatch.py`, which CLAUDE.md requires but this repository doesn't have, before any intensive job.
- [ ] Decide on the form of any write-up. Credit and relationship with Thong's paper need thought, since the draft is addressed to its author ("your draft"). Also decide whether to report the printed errors in Bulycheva 2018 to the author.

## Done

- [x] 2026-11-06: Read Thong 2026. Wrote the note, bibliography and synthesis.
- [x] 2026-11-06: Wrote `scripts/check_robin_identities.py`, testing the formulas of Thong and of the draft.
- [x] 2026-11-06: Verified the draft's citations on arXiv (bibliographic data only).
- [x] 2026-11-06: Showed the matrix boundary condition is self-adjoint with $\omega^2\ge0$; found the even-sector threshold $g_*$ and the homogeneous-mode relation (derivations 3).
- [x] 2026-11-06: Reviewed ten papers, in part where noted. Notes in `literature/notes/`.
- [x] 2026-11-06: **A1.** Compared [N2](2), (3), (9), (15), (16), (24), (26), (28) with the sources: all agree. Confirmed both claimed misprints in Bulycheva, and found more (`scripts/check_literature_inputs.py`, derivations 5).
- [x] 2026-11-06: Compared `kc_syk` with Maldacena–Stanford: agrees.
- [x] 2026-11-06: **B1.** Derived the relative cylinder from the super-Casimir (`scripts/check_supercasimir.py`, derivations 6).
- [x] 2026-11-06: **C2.** Resolved the $A$-branch modes (derivations 7).
