# Audit of the N=2 paper

- **Audited:** `research/pdf/n2_dirac_ladders.pdf`, sources in `research/tex/n2_dirac_ladders/`.
- **Version before the audit:** `research/pdf/n2_dirac_ladders_preaudit.pdf`, kept for the record.
- **Date:** 2026-11-06. There is no git repository yet, so this date stands in for a commit.
- **Script and log:** `scripts/audit_n2_paper.py`, `results/data/paper/audit_n2_paper.log`.

Section and equation numbers below refer to the corrected paper. Tags such as [B18] and [PSV] are as in `literature/notes/`.

## Verdict

**The main computational results hold.** The audit attacked them by routes the paper's own checks didn't use, and they survived:

- the reduced kernel as a Dirac correlator, by direct integration;
- the neutral kernel reducing to that same kernel, in position space;
- the closure of the supersymmetry algebra on the fermion;
- the breaking pattern $su(1,1|1)\to u(1|1)$;
- the reduced covariance [N2](12), now in all four entries at every root with $\omega<4$.

**Three statements in the pre-audit paper were wrong**, and are corrected:

1. "Supersymmetry is spontaneously broken below $g_*$" overstated the result. The Fock vacuum isn't invariant, but another state is.
2. "$\Phi_h$ has poles on $\chi>1$ only" is false. The conclusions drawn from it stand, by a different mechanism.
3. "[B18](4.99) is correct as printed" is false beyond the leading coefficient, which is all that had been tested.

There were also sign slips, a misattribution and several overstatements, listed below.

**What the paper doesn't establish** is unchanged. The super-centre integral of a superblock hasn't been done, so "one Dirac field replaces Thong's scalar" is still a conjecture for the full four-point function.

**A limit on this audit:** Claude wrote both the paper and the audit. It provides independent routes, not independent eyes, and shared conceptual errors could survive it. Nothing here has been externally reviewed.

## Errors found and corrected

### Wrong

**1. Spontaneous breaking below $g_*$** (abstract, Sec. 1 item 5, Sec. 7.2).

The paper said that below $g_*$ "supersymmetry is spontaneously broken", and the introduction called $2H-1\ge0$ a condition for a unitary multiplet. Working out the Fock space shows otherwise. For $\omega=H-\tfrac12<0$, the image $e^{-i\omega t}=e^{+i|\omega|t}$ has negative norm, so its coefficient in the field is a creation operator. Then $\bar Q=\sqrt{2|\omega|}\,af$, and
$$\tfrac12\{Q,\bar Q\}=|\omega|(N_a-N_f+1)=L_0-\tfrac12J_0,\qquad E_{\rm vac}=|\omega|=\tfrac12-H.$$

- The Fock vacuum, the lowest state of $L_0$, isn't invariant. That part was right, and the vacuum energy is confirmed as the constant in the algebra.
- The state with $N_a=0$, $N_f=1$ *is* invariant. It has $L_0=\tfrac12$, $J_0=1$, $L_0-\tfrac12J_0=0$, so it is a BPS state and the lowest state of the twisted Hamiltonian.
- The stationary quantization has positive norm at every $g$, so this isn't a unitarity bound.

The corrected Sec. 7.2 says that whether to call this spontaneous breaking depends on which Hamiltonian defines the vacuum. *Evidence:* part O. The matrices are tested. The operator content of $Q$ was read off from the norm and frequency, by hand.

**2. Poles of $\Phi_h$ "on $\chi>1$ only"** (Sec. 8.3, and `docs/derivations.md` 7.2).

This is false: $\mathrm{Res}_{h=1}\Phi^s_h=-2$ at $\chi=0.2,\ 0.5,\ 1.6$ (part I). The real difference between the regions is how the contour is closed. On $\chi>1$ one closes with $\Phi_h$. On $0<\chi<1$ one first uses $h\to1-h$ to pass to the block, following [MS], and the block has no pole there, so the zero of the measure removes the term. The conclusion stands: at $q=3$, $p=1$, $g=0.9999$, the family contributes $1.8\times10^{-6}$ at $\chi=0.5$ and $-0.26532$ at $\chi=1.6$, against the predicted $-0.26532$ (part L).

The wording "the four-point function does contain the mode" is also corrected. This is a term in the piecewise family decomposition, and doesn't imply that the full correlator jumps at $\chi=1$.

**3. "[B18](4.99) is correct as printed"** (App. A, and the Bulycheva note).

Only the leading coefficient had been tested. The printed factor $\sin\pi h-\sin2\pi\Delta/\cos\pi h$ should be $\sin2\pi\Delta-\sin\pi h/\cos\pi h$. It follows from the signs of (4.39) and (4.48), through (4.96). The two agree at $h=2\Delta$, so (4.100) stands, and they differ by $(-1)^n$ at $h=2\Delta+n$.

| $q=3$ | $\chi^{2\Delta}$ | corrected factor, full block | printed factor, full block | printed, bare $\chi^h$ |
|---|---|---|---|---|
| $\chi=0.2$ | 0.58480 | 0.58480 | 0.60214 | 0.59285 |
| $\chi=0.5$ | 0.79370 | 0.79370 | 0.86603 | 0.82237 |
| $\chi=0.8$ | 0.92832 | 0.92832 | 1.11915 | 0.98488 |

This is a new printed discrepancy. It doesn't affect the paper, which works from (4.106). [B18] says (4.99) was checked numerically for any $\chi$, so that check presumably used a correct form. *Evidence:* part P.

**4. Signs in Sec. 6.3.**

- "$\{q_A,\bar q_A\}=2\partial_\tau$, that is $2H_{\rm tw}$" had a sign slip. It is $-2H_{\rm tw}$ on $e^{-H_{\rm tw}\tau}$.
- "Anticommutator $2L_0-J_0$" depended on an unstated convention for $J_0$. It now states $J_0=+1$ on $\Psi$.
- The second line of (6.10) held only up to an unstated sign. Read literally from (6.4), $\delta\tilde{\bar a}=-(\partial_\tau A-\mathcal B)$. With $\tilde{\bar a}_F=-e^{-\tau/2}\bar a_F$, and both Grassmann parameters on the left, (6.10) is correct. With the other sign, $\{q_A,\bar q_A\}A=2\mathcal B$. This sign turns out to matter for the superpotential, below.

*Evidence:* parts E, E2.

**5. Citations.**

- The subtraction $F_{\rm sub}$ was attributed to [B18]. It comes from the earlier note, [N2](32).
- The Killing spinors were attributed to [CGS]. Their (A.7) is a different form, in other conventions, so the paper now says "compare".

### Overstated, now qualified

- **"Half of the supersymmetry."** It holds within one embedding. All four Killing spinors give odd symmetries of the boundary condition, with $\beta=+i$ for $A$ and $\beta=-i$ for $B$. What is lost is closure on AdS$_2$ isometries. The title still says "half-supersymmetric", which is a decision for you.
- **"$q_A$ is the residue of the CFT's total supertranslation."** Only the kernels match. Both $G_{-1/2}$ and $\bar G_{-1/2}$ have that kernel.
- **"These are Amsel and Marolf's superfields, written for $d=2$."** They work in $d=3,4,5$. Ours are the analogue.
- **"Every closed-form statement was tested."** Not so. App. B now lists what was done by hand only.
- **"Supported in one operator ordering"** (abstract). It now says a term in the family decomposition, supported on $\chi>1$.
- **"$Q/2$ should be $Q$"** in [B18]. The repair isn't unique, though $Q/2\to Q$ is the normalization of [FGMS].

Separately, the paper and notes used "her" for K. Bulycheva and "his" for Thong, which our sources don't establish. These are rewritten without pronouns.

## Gaps closed

| Statement | Before | After the audit |
|---|---|---|
| Closure on the fermion | Never tested | $\{q(\epsilon_Y),\bar q(\bar\epsilon_X)\}\Psi=-2\mathcal L_\xi\Psi+\alpha\Psi$ with $\alpha=+1,-1,0,0$ for $AA,BB,AB,BA$; $\{q,q\}=0$; closure on $\Phi_\psi$. 96 tests, part D |
| Fermion-to-boson maps, $\delta A$, $\delta\mathcal B$ | Only boson-to-fermion | Tested on Dirac modes of free frequency, part E2 |
| Neutral kernel on $(f,-\partial f)$ | Argued by hand through superspace | An exact component statement, Sec. 4.1. First row by direct quadrature, $\int G_bG_b=-\rho'$ from the chain formula, part N |
| Entries 12, 21, 22 of [N2](12) | Untested | Pass at every root with $\omega<4$, $q=3,5$, $g=0.9,1$, including exceptional roots, part K |
| Fermion part of the superpotential | "Inferred" | Derived, Sec. 5.3: $\delta(S_E+S_{\rm bdy})=\int d^2\theta(\Phi_+-c\Phi_-)\delta\Phi_-$, part N |
| $\mathcal N=1$ eigenvalue | By hand | From the [MSW] generators, part H |
| [B18](4.39) sign | Via (B.4), (B.5) | Direct integration of the definition (4.24), part G |
| $\omega^2\ge0$ | Argument only | $D_\psi$ real on the imaginary axis, and the eigenvalue relations, part N |

Two of these deserve comment.

**The superpotential is a real check of the structure.** The conditions (5.4) and (4.20) were found separately from the ladders. One superpotential now produces both, with the same $c$ and the same sign. The fermion terms agree only because of the sign in $\tilde{\bar a}_F$, which closure had fixed independently. A wrong structure would have given $b_F=-c\,a_F$.

**[PSV]'s mixed correlators need a normalization they don't state.** The matrix (4.25) acts on coefficients with no amplitudes extracted, so the normalized correlators of a column must carry $b_j^2$ of the $(3,4)$ pair. The literal reading of "and similarly for the others" below (4.5) gives $F_{21}/F_{12}=\delta^2$, which violates the symmetry of the covariance. The column normalization gives 1. This is an inference from that symmetry, and without it entries 21 and 12 of [N2](12) would fail.

## Claim by claim

**Status:** *confirmed* means tested by a route other than the derivation, such as direct integrals, independent mode functions, or the printed source. *Algebra* means a formula against a formula, or exact algebra on test values. *Hand* means derived only. All numerical statements hold at the tested parameter values, and aren't proofs.

| Claim | Sec. | Status | Evidence |
|---|---|---|---|
| Kernels $k_{A,S}$, spectra, multiplets | 2 | Confirmed | Transcribed from [PSV], [B18], [FGMS]; `check_literature_inputs.py` |
| Generators and Casimir repairs | 3 | Confirmed | Closure and centrality on a Grassmann algebra; reproduces (A.10)–(A.12), (A.16); `check_supercasimir.py`, part H |
| $WCW^{-1}=\partial_\tau^2$, no odd relative coordinate | 3 | Confirmed | Same; counting argument by hand |
| $h(h-d_s)$ for $\mathcal N=0,1,2$ | 3 | Confirmed | Part H for $\mathcal N=1$ |
| Reduced kernel, eigenvalue $k_p(H)$ | 4 | Confirmed | Direct $y$-integration, both signs of $x$, $q=3,5$, part B |
| Seed and kernel as Dirac boundary correlators | 4 | Confirmed | Free spectral sum against closed forms, part C; (4.10) by quadrature |
| Neutral kernel equals $K^{(0)}$ on $(f,-\partial f)$ | 4 | Confirmed | Part N. The superspace explanation of it is hand only |
| $\hat{\mathsf k}_p=\mathsf D_p\mathsf M_p$; Ward identity | 5 | Confirmed | Entrywise against [PSV]; real, complex and imaginary $\omega$; $q=3,5,7$ |
| Boundary maps, $\Phi_+=c\Phi_-$ | 5 | Confirmed | Modes, `check_bulk_susy_modes.py`; algebra in Prop. 5.1 |
| Flux and superpotential | 5 | Confirmed in part | Wronskians, Dirac bilinear, Grassmann algebra tested. Boundary term $\bar a_Fb_F$ and Euclidean continuation by hand |
| Killing spinors, intertwiners, two embeddings | 6 | Confirmed | Symbolic, `check_bulk_susy_algebra.py` |
| Closure on $\Psi$ and $\Phi_\psi$ | 6 | Confirmed | Part D |
| Boundary action of $q_A,\bar q_A$ | 6 | Confirmed | Parts E, E2; closes in matrix form |
| $su(1,1\vert1)\to u(1\vert1)$ for the given domain | 6 | Confirmed | Residuals on modes; $E_+$ on bulk modes and its mismatch, part F |
| Mixed anticommutator $-4\omega$; four charges don't close | 6 | Algebra | On modes. The algebra they generate isn't classified |
| $\omega^2\ge0$ | 7 | Confirmed | Part N, with self-adjointness |
| Threshold $g_*$ and its values | 7 | Algebra | Root-finding at $q=3,5,7$. That only this root reaches $H=\tfrac12$ is numerical |
| $E_{\rm vac}$, non-invariant vacuum, BPS state | 7 | Algebra | Part O; operator content by hand; sum per multiplet, no regularization of the full series |
| $F_g(z)$ smooth across $g_*$ | 7 | Confirmed | Contour against root sum at $g=0.15$, $q=3$, `check_reconstruction.py` |
| Continued against stationary covariance | 7 | Algebra | $q=3,5$ |
| Reconstruction (8.3), (8.5) | 8 | Confirmed | Numerically at finite $g$; residues against [PSV](4.26) at $g=0.7,1$; seed gives $z^\delta$ |
| Reduced covariance (8.6), entry 11 | 8 | Confirmed | From [PSV] and [T](C.8), every root $\omega<4$, $g=0.5,0.9,0.99,1$ ($p=1$), $g=0.9,0.99,1$ ($p=0$) |
| Reduced covariance, other entries | 8 | Confirmed, given the normalization | Part K |
| $\mathcal I_a(h)$ for the even sector | 8 | Quadrature | With an analogy to Thong's weights, part M. No analytic derivation |
| Odd seed reduces correctly | 8 | Confirmed | Direct centre integral, part J. Even part has no such test |
| Exceptional families, step-function lemma | 8 | Confirmed | Lemma proved and tested; limit in part L; entries at $h=1,2$ in part K |
| Finite part of the $S$ multiplet | 8 | Algebra | Coefficients only. Not the physical answer, as the paper says |
| Complex SYK relation | 9 | Algebra | $q=4,6,8,12$, complex $h$. Seed and blocks not reduced |
| Supersymmetric line, $g=1-r$ | 7 | Hand | A proposal, labelled so |
| Printed discrepancies in [B18] | A | Confirmed | Against the PDF, by eye, for (4.39), (4.48), (4.49), (4.88), (4.90), (4.92), (4.96), (4.98), (4.99), (A.3), (A.5), (A.7), (A.9), (A.15); by computation for all |

## Not established

1. **The super-centre integral of a superblock.** This remains the main gap.
2. **Reduction of the summed series.** The direct test is only at $g=0$, in the odd sector, for $q=3,5$. Convergence isn't proved.
3. **The even-sector reduction.** $\mathcal I_a$ comes from quadrature, and the seed's centre integral diverges, so there is no direct test.
4. **A direct two-rung integral at order $g^2$.**
5. **The Fock charges from a mode expansion.** Their form is read from the norm and frequency. The role of the BPS state in the CFT is unknown.
6. **The boundary term and continuation in the flux,** and the superspace reduction of [B18](4.42). Both by hand.
7. **The supersymmetric line.** No component reduction of $S_0$, no reflection positivity, nothing at order $1/N$.
8. **The algebra of the four odd symmetries,** and whether a massive Type II condition exists.
9. **Generality in $q$.** Almost everything is at $q=3,5$ or $3,5,7$.
10. **Novelty.** The literature search isn't exhaustive. The Amsel–Marolf PDF isn't in `papers/`, and the Maldacena–Stanford PDF is v1.
11. **Transcription.** Formulas were taken from PDF text extraction. Only the [B18] equations in the table above were also checked by eye.

## How much the passing tests mean

**A passing count doesn't measure truth.** The pre-audit paper passed all 1124 of its checks and still contained the errors above, because those statements weren't what the tests tested. The pole statement had no test on $0<\chi<1$. (4.99) was tested at one coefficient. The paper's own closure check never touched the fermion. Most earlier tests compare one closed form with another, written by the same author.

**The audit script made mistakes too.**

- Its first version had a "BPS check" that compared a formula with itself. It passed, and meant nothing. Part O replaces it, and it was while replacing it that error 1 came to light.
- It produced false failures from numerics: slow power-law tails in parts A and B, cancellation noise in part G, and a $u^{-0.8}$ endpoint in part N. These were fixed by variable changes, hand-folded principal values and higher precision. In part N the diagnosis was confirmed against exact Beta functions.
- Tolerances are per test, and some are coarse on purpose. The boundary fit in part A agrees to about 1 part in $10^5$, because it truncates an expansion.
- The 34 symmetry checks in part K, $F_{21}/F_{12}$, count in the total but print as `[INFO]` lines rather than `[PASS]` lines.

**Independence is partial.** The script has its own kernels, mode functions and direct integrals, but imports the [PSV] formulas from `check_literature_inputs.py` and the Grassmann helpers from `check_supercasimir.py`.

## Reproducing

```bash
python3 scripts/audit_n2_paper.py --dps 20               # all parts
python3 scripts/audit_n2_paper.py --dps 20 --parts KNOP  # selected parts; lowercase e is E2
```

- **Result:** 472/472 checks passing: 438 `[PASS]` lines, plus the 34 symmetry checks of part K.
- **Earlier checks, rerun afterwards:** 1124/1124 passing, unchanged. Logs in `results/data/paper/rerun_after_audit/`. The earlier logs in `results/data/paper/` are left as they were.
- **Parameters:** 20 working digits, higher where the script says so. $q=3,5$ in most parts, and $q=7$ for the boson rail. $g=0.9,1$ in part K. Tolerances are per test, in the script. Part K uses $10^{-6}$, part J $10^{-5}$, and part B $10^{-8}$.
- **Cost:** a few minutes on one core, with scalar `mpmath` arithmetic and no large arrays. It isn't an intensive job in the sense of `CLAUDE.md`, so it didn't run under the memory watchdog, which isn't in the repository yet. No random numbers are used.
- **Rebuilding the paper:** in `research/tex/n2_dirac_ladders/`, run `pdflatex main.tex` three times, and copy `main.pdf` to `research/pdf/n2_dirac_ladders.pdf`.

The living documents were updated to match: `docs/derivations.md` (Secs. 7.2, 7.3, 10.3–10.6, 11), `docs/research_questions.md`, `docs/todo.md`, `literature/synthesis.md`, and the Bulycheva, PSV and Amsel–Marolf notes.

## Addendum, 2026-11-07

The first item under "Not established", the super-centre integral of a superblock, has since been done. See paper Sec. 8.3, `docs/derivations.md` Sec. 12, and `scripts/check_super_centre.py` (182/182). Each superblock reduces to one exponential, and the four-point function reduces to the Dirac propagator. The audit text above is left as it was.
