# Research questions

**Project goal:** extend the Robin-ladder construction of Thong (arXiv:2610.31989) to N=2 supersymmetric SYK, starting from the AI-generated draft `research/pdf/n2_susy_generalization_note.pdf`, and find out which of its claims hold.

Evidence and derivations are in [`derivations.md`](derivations.md). The order of work is in [`todo.md`](todo.md).

*Last updated 2026-11-06, after the first literature review.*

**Summary**

| | Question | Status |
|---|---|---|
| A1 | Physical input against the literature | **Answered: agrees** |
| A2 | Reconstruction at finite $g$ | Partly answered |
| A3 | Independent data | Partly answered |
| B1 | Is there a super-cylinder? | **Answered: bosonic cylinder, $d_s=0$** |
| B2 | Centre reduction of superblocks | Component level done; super level open |
| B3 | Mixed continuation | Open |
| C1 | Is $g=1$ the stationary covariance? | Mostly answered |
| C2 | $A$-branch $H=1$ modes | **Answered: step functions on $\chi>1$** |
| C3 | Minimal field content | Conjecture: one Dirac field |
| C4 | Surviving symmetry | Open, with a concrete next test |
| C5 | Physical $H=1$ coefficient | Open |
| D1 | Simpler cases | Answered for complex SYK |
| D2 | Is $g\neq1$ meaningful? | Proposal: a supersymmetric line |
| D3, D4 | Charge; beyond four points | Open |
| D5 | What's new? | Open; search not done |

## A. Is the draft right?

**A1. Does the physical input match the literature?** Do [N2](2), (3), (15) and the overlaps behind (26) really agree with Peng–Spradlin–Volovich and Bulycheva?

*Status: answered, yes (derivations 5.1, 5.3). The draft's two claimed misprints in Bulycheva are both confirmed (5.2), and we found further ones.*

**A2. Is the four-point reconstruction right at finite $g$?** The representation [N2](27) reproduces orders $g^0$ and $g^1$. Order $g^2$ against a direct two-rung integral would be a real test, because it probes the kernel structure beyond the seed and one rung.

*Status: partly answered. The residues (28) agree with [PSV] at finite $g$ (derivations 5.4). The contour form and the $g^2$ integral are still to do.*

**A3. Does the construction reproduce independent data?** Compare with the large-$N$ four-point results of the source papers, or with finite-$N$ numerics if a clean observable exists.

*Status: partly answered. [N2](12) agrees with [PSV]+[T](C.8) in the 11 entry, for $\omega<4$ (derivations 5.4). Other entries, and any finite-$N$ test, remain.*

## B. Can Thong's derivation really be carried over?

**B1. Is there a super-cylinder?** What does zero total supermomentum do to the $su(1,1|1)$ pair Casimir? Does this derive the relative cylinder, including $\omega=H$ and shadow $-H$, rather than inferring it from kernels?

*Status: answered (derivations 6). The relative space is purely bosonic, $f=f(\langle12\rangle)$. $C\to L_0^2$, so $WCW^{-1}=\partial_\tau^2$ with $a=\delta$. In general $d_s=1-\mathcal N/2$. Not a "super-cylinder" in the sense of having odd directions.*

**B2. What is the N=2 centre reduction of a superblock?** What are the analogues of [T](C.5)–(C.8), across all orderings? Do the reduced families stay nonzero and independent?

*Status: partly answered. Ordinary centre integrals of component blocks: $\mathcal I_a(h)=\pi\cot(\pi h/2)/[2\omega\mathcal N]$ found by quadrature, plus exact step-function families (derivations 7). **Open, and now the main gap:** the super-centre integral of a superblock, and whether $\mathcal T_2$ is its inverse.*

**B3. Is the mixed $\psi$–$b$ continuation consistent?** Do the prescriptions for the IR-divergent fermion rail and UV-divergent boson rail allow the centre shifts and integral exchanges Thong requires?

*Status: open. The super formulation might avoid this, as one rail $G(\langle12\rangle)$ replaces two.*

## C. What does the auxiliary system mean?

**C1. Is the reduced $g=1$ ladder the stationary covariance?** The matrix boundary problem is self-adjoint with $\omega^2\ge0$ (derived). The even-parity seed fails $a+J>0$ by one homogeneous mode, and there is a zero mode at a formal $g_*<1$.

*Status: mostly answered. Mode by mode, yes in the 11 entry for $\omega<4$, without needing continuity through $g_*$ (derivations 3.5, 5.4). The summed statement and other entries remain.*

**C2. What are the $A$-branch $H=1$ modes?** They carry nonzero weight in the auxiliary covariance, but zero weight in the draft's chiral–antichiral reconstruction.

*Status: answered (derivations 7). The block coefficient vanishes while $\mathcal I$ has a pole. The family is a finite step function supported only on $\chi>1$, which the draft's $0<z<1$ can't see, and its centre integral is exactly the mode. [T](3.13) holds. This matches [FGMS]'s expectation that this multiplet gives IR operators, not the super-Schwarzian. **New question:** what does a contribution in only one ordering mean physically?*

**C3. What is the minimal field content?** Two scalars for each parity, and a Dirac field for the charged channels, are linked by first-order radial equations. What is the N=2 replacement for Thong's "single field"?

*Status: conjecture (derivations 6). One Dirac field with $b_F=g\sigma L_qa_F$. The matrix is what ordinary, rather than super, centre integration gives: it keeps both bottom and top components of each multiplet. Needs B2 to prove.*

**C4. Which symmetry survives?** Thong's Robin condition leaves $\mathbb R\times SO(d)$. The draft shows that the N=2 boundary domain does not preserve the full AdS isometry algebra [N2](35), so not all four AdS N=2 supercharges. Does a smaller subalgebra survive?

*Status: open, with a concrete test. [CGS] shows that a plain Robin condition can keep 2 of 4 supercharges (massless case), and that the supercharges have global energy $\pm\tfrac12$. Test the draft's domain against the Killing spinor directly, with $m_F=\mu$, instead of only $E_+$. [CGS] Type II, with a derivative condition on the scalar and an algebraic one on the fermion, looks like the structure of [N2](20), (23), but [CGS] covers only $m=0$.*

**C5. What is the physical content of the $H=1$ multiplet?** Connect the finite part [N2](33) to the super-Schwarzian of Fu et al. The draft states that the $\varphi_1$ coefficient depends on the subtraction.

*Status: open. Now sharper: this concerns the $S$ multiplet only (derivations 7.3). Read [MS] Sec. 3.3 and [FGMS] (5.27)–(5.38). [PSV] and [B18] both left it undone, so this may be new.*

## D. Scope questions

**D1. Would simpler cases help?**

*Status: answered for complex SYK. Its symmetric channel is exactly Thong's $J=0$ sector, $k^S\mathcal R_0=4^{-\nu}$ (derivations 9). A single scalar with a parity-dependent $\kappa$, $a+J<0$ for $q>4$, and an $h=1$ mode. It's the best first test for going past Thong's validity range. N=1 has an odd relative coordinate, so it tests something else, and isn't simpler.*

**D2. Is $g\neq1$ meaningful?** Is there an N=2 analogue of the Gross–Rosenhaus line preserving supersymmetry, or is $g$ purely a formal device?

*Status: proposal (derivations 8). A superspace bilocal kinetic term gives $g=1-b/\bar b$, and scales every kernel entry uniformly at leading large $N$, which is exactly the draft's $g$. **New question:** is $g_*$ then a physical transition on the line, where the lowest even-sector operator reaches $h=\tfrac12$? It has no $\mathcal N=0$ analogue in [GR].*

**D3. What about nonzero charge?** The draft works only at the charge-symmetric saddle. Does the construction survive with a chemical potential?

*Status: open; not a priority yet. [GKST] App. B gives a worked example of $\theta$-dependent channel mixing with protected $h=1,2$.*

**D4. What about beyond four points?** Thong raises cubic couplings and higher-point data as open even for $\mathcal N=0$. The draft notes the six-point case is not done.

*Status: open; later.*

**D5. What is new relative to the literature?** Before calling anything new, search for prior or parallel work. The draft states "no priority claim is made" and has not done this.

*Status: open, and now more pressing, given Secs. 6–8 of the derivations. **None of those results has been checked against the literature beyond these ten papers.***
