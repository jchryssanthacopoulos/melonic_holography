# Research questions

**Project goal:** extend the Robin-ladder construction of Thong (arXiv:2610.31989) to N=2 supersymmetric SYK, starting from the AI-generated draft `research/pdf/n2_susy_generalization_note.pdf`, and find out which of its claims hold.

Evidence and derivations are in [`derivations.md`](derivations.md). The order of work is in [`todo.md`](todo.md).

*Last updated 2026-11-07, after the super-centre integral (derivations 12). Earlier: the literature review, the paper draft (`research/pdf/n2_dirac_ladders.pdf`) and its audit (`research/notes/n2_dirac_ladders_audit.md`).*

**Summary**

| | Question | Status |
|---|---|---|
| A1 | Physical input against the literature | **Answered: agrees** |
| A2 | Reconstruction at finite $g$ | Partly answered |
| A3 | Independent data | Partly answered |
| B1 | Is there a super-cylinder? | **Answered: bosonic cylinder, $d_s=0$** |
| B2 | Centre reduction of superblocks | **Answered: one exponential per superblock; $F_g\to$ Dirac propagator** |
| B3 | Mixed continuation | Open |
| C1 | Is $g=1$ the stationary covariance? | Mostly answered |
| C2 | $A$-branch $H=1$ modes | **Answered: step functions on $\chi>1$** |
| C3 | Minimal field content | **Answered at four points: one Dirac field; scalars are components** |
| C4 | Surviving symmetry | **Answered: $su(1,1|1)\to u(1|1)$, $4|4\to2|2$** |
| C5 | Physical $H=1$ coefficient | Open |
| D1 | Simpler cases | Answered for complex SYK |
| D2 | Is $g\neq1$ meaningful? | Proposal: a supersymmetric line; $g_*$ is a threshold for the Fock vacuum, not a CFT transition |
| D3, D4 | Charge; beyond four points | Open |
| D5 | What's new? | Partial search: nothing found on the N=2 extension. Closest precedent: Amsel–Marolf |

## A. Is the draft right?

**A1. Does the physical input match the literature?** Do [N2](2), (3), (15) and the overlaps behind (26) really agree with Peng–Spradlin–Volovich and Bulycheva?

*Status: answered, yes (derivations 5.1, 5.3). The draft's two claimed misprints in Bulycheva are both confirmed (5.2), and we found further ones.*

**A2. Is the four-point reconstruction right at finite $g$?** The representation [N2](27) reproduces orders $g^0$ and $g^1$. Order $g^2$ against a direct two-rung integral would be a real test, because it probes the kernel structure beyond the seed and one rung.

*Status: partly answered. The residues (28) agree with [PSV] at finite $g$ (derivations 5.4). The contour form and the $g^2$ integral are still to do.*

**A3. Does the construction reproduce independent data?** Compare with the large-$N$ four-point results of the source papers, or with finite-$N$ numerics if a clean observable exists.

*Status: partly answered. [N2](12) agrees with [PSV]+[T](C.8) in all four entries, at every root with $\omega<4$, for $q=3,5$ and $g=0.9,1$, including the exceptional roots, once [PSV]'s mixed correlators are normalized by the $(3,4)$ pair (audit part K; derivations 5.4, 11). Any finite-$N$ test remains.*

## B. Can Thong's derivation really be carried over?

**B1. Is there a super-cylinder?** What does zero total supermomentum do to the $su(1,1|1)$ pair Casimir? Does this derive the relative cylinder, including $\omega=H$ and shadow $-H$, rather than inferring it from kernels?

*Status: answered (derivations 6). The relative space is purely bosonic, $f=f(\langle12\rangle)$. $C\to L_0^2$, so $WCW^{-1}=\partial_\tau^2$ with $a=\delta$. In general $d_s=1-\mathcal N/2$. Not a "super-cylinder" in the sense of having odd directions.*

**B2. What is the N=2 centre reduction of a superblock?** What are the analogues of [T](C.5)–(C.8), across all orderings? Do the reduced families stay nonzero and independent?

*Status: answered (derivations 12, paper Sec. 8.3). The super-centre integral of a superblock is a single one-sided exponential $e^{-H|t|}$, with weight $-4\pi\tan(\pi H/2)$ ($A$) or $4\pi\cot(\pi H/2)$ ($S$), and the four-point function reduces to the Dirac propagator, $8\pi\mathcal N_\delta L_qb_\psi^2\,\mathfrak G_{F,g}\Pi$. The reconstruction built on $\mathcal T_2$ is its inverse. Ordinary centre integrals of component blocks: $\mathcal I_a(h)=\pi\cot(\pi h/2)/[(2h-1)\mathcal N]$ found by quadrature, plus exact step-function families (derivations 7). **Still open:** term-by-term reduction of the infinite sum, and a superspace characterization of the continued superfamilies.*

**B3. Is the mixed $\psi$–$b$ continuation consistent?** Do the prescriptions for the IR-divergent fermion rail and UV-divergent boson rail allow the centre shifts and integral exchanges Thong requires?

*Status: open. The super formulation might avoid this, as one rail $G(\langle12\rangle)$ replaces two.*

## C. What does the auxiliary system mean?

**C1. Is the reduced $g=1$ ladder the stationary covariance?** The matrix boundary problem is self-adjoint with $\omega^2\ge0$ (derived). The even-parity seed fails $a+J>0$ by one homogeneous mode, and there is a zero mode at a formal $g_*<1$.

*Status: mostly answered. Mode by mode, yes in all four entries for $\omega<4$ at the tested $q,g$, without needing continuity through $g_*$ (derivations 3.5, 5.4, 11). The summed statement remains; it has a direct test only at $g=0$ in the odd sector.*

**C2. What are the $A$-branch $H=1$ modes?** They carry nonzero weight in the auxiliary covariance, but zero weight in the draft's chiral–antichiral reconstruction.

*Status: answered (derivations 7). The block coefficient vanishes while $\mathcal I$ has a pole. In the piecewise decomposition, the family is a finite step function supported only on $\chi>1$, which the draft's $0<z<1$ can't see, and its centre integral is exactly the mode. [T](3.13) holds. This is a term in a decomposition, not a jump of the full correlator. This matches [FGMS]'s expectation that this multiplet gives IR operators, not the super-Schwarzian. **New question:** what does a contribution in only one ordering mean physically?*

**C3. What is the minimal field content?** Two scalars for each parity, and a Dirac field for the charged channels, are linked by first-order radial equations. What is the N=2 replacement for Thong's "single field"?

*Status: answered at the level of four-point data (derivations 10, 12). The reduced ladder is one Dirac propagator with a boundary mass. The neutral $2\times2$ correlator is that Dirac correlator with kinematic projectors, at all complex $\omega$. The scalars complete the Dirac field to a long AdS$_2$ multiplet with $H=\mp\mu$. And the super-centre integral of the full neutral four-point function is the same Dirac propagator, family by family. As in [T], convergence of the family sums isn't proved.*

**C4. Which symmetry survives?** Thong's Robin condition leaves $\mathbb R\times SO(d)$. The draft shows that the N=2 boundary domain does not preserve the full AdS isometry algebra [N2](35), so not all four AdS N=2 supercharges. Does a smaller subalgebra survive?

*Status: answered for the given boundary domain (derivations 10.4, paper Sec. 6). $su(1,1|1)\to u(1|1)$ in each embedding: besides $E_\pm$, two of the four Killing supercharges break, and the two that survive close on $L_0-\tfrac12J_0$. The other two still give odd symmetries of the boundary condition in the other embedding, so "half" refers to any one $su(1,1|1)$. **Still open:** whether a massive analogue of [CGS] Type II could preserve more, and the algebra of the four odd symmetries from the two embeddings.*

**C5. What is the physical content of the $H=1$ multiplet?** Connect the finite part [N2](33) to the super-Schwarzian of Fu et al. The draft states that the $\varphi_1$ coefficient depends on the subtraction.

*Status: open. Now sharper: this concerns the $S$ multiplet only (derivations 7.3). Read [MS] Sec. 3.3 and [FGMS] (5.27)–(5.38). [PSV] and [B18] both left it undone, so this may be new.*

## D. Scope questions

**D1. Would simpler cases help?**

*Status: answered for complex SYK. Its symmetric channel is exactly Thong's $J=0$ sector, $k^S\mathcal R_0=4^{-\nu}$ (derivations 9). A single scalar with a parity-dependent $\kappa$, $a+J<0$ for $q>4$, and an $h=1$ mode. It's the best first test for going past Thong's validity range. N=1 has an odd relative coordinate, so it tests something else, and isn't simpler.*

**D2. Is $g\neq1$ meaningful?** Is there an N=2 analogue of the Gross–Rosenhaus line preserving supersymmetry, or is $g$ purely a formal device?

*Status: proposal (derivations 8). A superspace bilocal kinetic term gives $g=1-b/\bar b$, and scales every kernel entry uniformly at leading large $N$, which is exactly the draft's $g$. $g_*$ is where the Fock vacuum stops being invariant under the surviving charges, and the boundary of validity of the reduced bosonic description. It is not a transition of $F_g(z)$ on $0<z<1$. Below it, a state with one Dirac mode filled is invariant, with $L_0=\tfrac12$, $J_0=1$ (derivations 10.5). The super-centre reduction to the Dirac propagator has no boundary at $g_*$ (derivations 12). **New question:** does that BPS state have a role in the CFT?*

**D3. What about nonzero charge?** The draft works only at the charge-symmetric saddle. Does the construction survive with a chemical potential?

*Status: open; not a priority yet. [GKST] App. B gives a worked example of $\theta$-dependent channel mixing with protected $h=1,2$.*

**D4. What about beyond four points?** Thong raises cubic couplings and higher-point data as open even for $\mathcal N=0$. The draft notes the six-point case is not done.

*Status: open; later.*

**D5. What is new relative to the literature?** Before calling anything new, search for prior or parallel work. The draft states "no priority claim is made" and has not done this.

*Status: partial. Web searches (2026-11-06) found no prior N=2 extension of Thong's construction, but they weren't exhaustive. The boundary-superfield structure and the half-supersymmetry pattern have a precedent in Amsel–Marolf for $d=3,4,5$ (see its note). The rest of derivations 6–10 hasn't been checked against the literature beyond the papers we read.*
