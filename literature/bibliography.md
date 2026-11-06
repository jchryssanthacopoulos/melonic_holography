# Bibliography

**Status** says how far each entry has been read:

- **read**: the PDF is in `papers/` and has a note in `literature/notes/`.
- **part**: the same, but read only in part. The note says which sections, and claims should stay within them.
- **background**: cited for context only, not in `papers/`, and not read. Claims attributed to it are second-hand.

**Metadata source** is one of: *arXiv* (checked on the arXiv abstract page, 2026-11-06), *PDF* (arXiv number read from the PDF in `papers/`, 2026-11-06), *memory* (the journal reference, flagged for checking) or *Thong refs* (copied from the reference list of Thong 2026, not independently checked).

Bracketed numbers such as `[T36]` are the reference numbers in Thong 2026.

## Core paper

| Key | Citation | Status | Metadata source |
|---|---|---|---|
| `thong_2026` | K. H. Thong, *Melonic ladder holography with a single free bulk field in AdS<sub>d+1</sub>*, arXiv:2610.31989v1 [hep-th] (22 Oct 2026), prepared for JHEP. PDF: `papers/thong_2026.pdf`. Note: [`notes/thong_2026.md`](notes/thong_2026.md). | read | arXiv |

## N=2 SYK and its sources

| Key | Citation | Why it matters | Status | Metadata source |
|---|---|---|---|---|
| `peng_spradlin_volovich_2017` | C. Peng, M. Spradlin, A. Volovich, *Correlators in the N=2 Supersymmetric SYK Model*, JHEP 10 (2017) 202, arXiv:1706.06078v2. [T36] [Note](notes/peng_spradlin_volovich_2017.md) | Component kernels, spectra, four-point functions. The draft agrees. | read | arXiv, PDF |
| `bulycheva_2018` | K. Bulycheva, *N=2 SYK model in the superspace formalism*, JHEP 04 (2018) 036, arXiv:1801.09006v2. [Note](notes/bulycheva_2018.md) | Superblocks, Casimir, measure. Several printed errors, listed in the note. Gives B1. | part (Secs. 1–4, 7, Apps. A–G) | arXiv, PDF |
| `fu_gaiotto_maldacena_sachdev_2016` | W. Fu, D. Gaiotto, J. Maldacena, S. Sachdev, *Supersymmetric SYK models*, Phys. Rev. D 95 (2017) 026009, arXiv:1610.08917. [Note](notes/fu_gaiotto_maldacena_sachdev_2016.md) | Defines the model; super-Schwarzian; the two $H=1$ multiplets. | part (Secs. V–VII) | arXiv, PDF |
| `correa_giraldo-rivera_silva_2019` | D. H. Correa, V. I. Giraldo-Rivera, G. A. Silva, *Supersymmetric mixed boundary conditions in AdS<sub>2</sub> and DCFT<sub>1</sub> marginal deformations*, JHEP 03 (2020) 010, arXiv:1910.04225v3. [Note](notes/correa_giraldo-rivera_silva_2019.md) | Supersymmetric boundary conditions, massless case. Bears on C4. | part (Secs. 1, 2, 3.1, 3.2, 5, App. A) | arXiv, PDF |
| `bulycheva_2017` | K. Bulycheva, *A note on the SYK model with complex fermions*, JHEP 12 (2017) 069, arXiv:1706.07411v1. [T35] [Note](notes/bulycheva_2017.md) | Symmetric channel; $U(1)$ mode. A different paper from `bulycheva_2018`. | part (Sec. 1, parts of 4, 7) | Thong refs, PDF |
| `murugan_stanford_witten_2017` | J. Murugan, D. Stanford, E. Witten, *More on supersymmetric and 2d analogs of the SYK model*, JHEP 08 (2017) 146, arXiv:1706.05362v1. [Note](notes/murugan_stanford_witten_2017.md) | N=1 Casimir and kernels. | part (targeted; Sec. 2, start of 3.6) | PDF; journal from memory |
| `gu_kitaev_sachdev_tarnopolsky_2019` | Y. Gu, A. Kitaev, S. Sachdev, G. Tarnopolsky, *Notes on the complex Sachdev-Ye-Kitaev model*, JHEP 02 (2020) 157, arXiv:1910.14099v3. [Note](notes/gu_kitaev_sachdev_tarnopolsky_2019.md) | Charged spectrum; protected modes; Dirac bulk for the entropy. | part (targeted; Sec. 5.1, App. B) | PDF; journal from memory |

## SYK inputs that Thong's derivation uses directly

| Key | Citation | Used for | Status | Metadata source |
|---|---|---|---|---|
| `maldacena_stanford_2016` | J. Maldacena, D. Stanford, *Remarks on the Sachdev-Ye-Kitaev model*, Phys. Rev. D 94 (2016) 106002, arXiv:1604.07818. [T1] PDF: `papers/maldacena_standford_2016.pdf`, which is **v1, titled "Comments on..."**, with numbering that differs from what Thong cites. [Note](notes/maldacena_stanford_2016.md) | Conformal kernel, OPE residues, h=2 sector. | part (Secs. 3.2.1–3.2.5) | Thong refs, PDF |
| `gross_rosenhaus_2017` | D. J. Gross, V. Rosenhaus, *A line of CFTs: from generalized free fields to SYK*, JHEP 07 (2017) 086, arXiv:1706.07015v2. [T17] [Note](notes/gross_rosenhaus_2017.md) | The line 0 ≤ g < 1 and its OPE residues. Basis for D2. | part (Secs. 1, 2, 3.1) | Thong refs, PDF |
| `giombi_klebanov_tarnopolsky_2017` | S. Giombi, I. R. Klebanov, G. Tarnopolsky, *Bosonic tensor models at large N and small ε*, Phys. Rev. D 96 (2017) 106014, arXiv:1707.03866. [T18] | Chain integral and kernel eigenvalue for d ≥ 2. | background | Thong refs |
| `benedetti_delporte_2021` | D. Benedetti, N. Delporte, *Remarks on a melonic field theory with cubic interaction*, JHEP 04 (2021) 197, arXiv:2012.12238. [T24] | AR model conventions. Thong uses half of the eigenvalue printed in its Eq. (6.19). | background | Thong refs |
| `ishibashi_wald_2004` | A. Ishibashi, R. M. Wald, *Dynamics in nonglobally hyperbolic static space-times. 3. Anti-de Sitter space-time*, Class. Quant. Grav. 21 (2004) 2981, hep-th/0402184v4. [T64] [Note](notes/ishibashi_wald_2004.md) | Self-adjoint extensions and the positivity theorem, single-field. | part (Sec. 3.3) | Thong refs, PDF |

## Alternative single-field and bulk constructions for SYK

| Key | Citation | Relevance | Status | Metadata source |
|---|---|---|---|---|
| `das_ghosh_jevicki_suzuki_2018` | S. R. Das, A. Ghosh, A. Jevicki, K. Suzuki, *Three Dimensional View of Arbitrary q SYK models*, JHEP 02 (2018) 162, arXiv:1711.09839. [T9] | The "DGJS" interval problem, mapped to Thong's radial problem in appendix E. | background | Thong refs |
| `das_jevicki_suzuki_2017` | S. R. Das, A. Jevicki, K. Suzuki, *Three Dimensional View of the SYK/AdS Duality*, JHEP 09 (2017) 017, arXiv:1704.07208. [T10] | As above. | background | Thong refs |
| `vegh_2024` | D. Vegh, *Quantizing the folded string in AdS<sub>2</sub>*, arXiv:2409.06663. [T12] | Folded string. Source of the coordinate map in appendix F. | background | Thong refs |
| `vegh_2025` | D. Vegh, *A folded string dual for the Sachdev-Ye-Kitaev model*, arXiv:2509.05435. [T13] | Folded-string dual; Thong's stated inspiration. | background | Thong refs |
| `thong_vegh_2026` | K. H. Thong, D. Vegh, *Quantum mechanical bootstrap without inequalities: SYK bilinear spectrum*, arXiv:2604.26007. [T14] | Fractional-power Robin conditions on an interval. | background | Thong refs |
| `gross_rosenhaus_2017_cubic` | D. J. Gross, V. Rosenhaus, *The Bulk Dual of SYK: Cubic Couplings*, JHEP 05 (2017) 092, arXiv:1702.08016. [T8] | Cubic couplings: the natural next test beyond four points. | background | Thong refs |

## Robin boundary conditions and double-trace deformations

| Key | Citation | Relevance | Status | Metadata source |
|---|---|---|---|---|
| `witten_2001` | E. Witten, *Multitrace operators, boundary conditions, and AdS/CFT correspondence*, hep-th/0112258. [T30] | Quadratic boundary terms as mixed boundary conditions. | background | Thong refs |
| `gubser_klebanov_2003` | S. S. Gubser, I. R. Klebanov, *A Universal result on central charges in the presence of double trace deformations*, Nucl. Phys. B 656 (2003) 23, hep-th/0212138. [T31] | Double-trace resolvent. | background | Thong refs |
| `hartman_rastelli_2008` | T. Hartman, L. Rastelli, *Double-trace deformations, mixed boundary conditions and functional determinants in AdS/CFT*, JHEP 01 (2008) 019, hep-th/0602106. [T32] | Double-trace resolvent. | background | Thong refs |
| `morley_taylor_winstanley_2021` | T. Morley, P. Taylor, E. Winstanley, *Quantum field theory on global anti-de Sitter space-time with Robin boundary conditions*, Class. Quant. Grav. 38 (2021) 035009, arXiv:2004.02704. [T33] | Boundary-chain expansions and a boundary-condition instability in global AdS<sub>4</sub>. | background | Thong refs |
| `wang_yang_2026` | Y. Wang, J. Yang, *Robin holography in AdS and BTZ: double-trace RG flow and exceptional points*, JHEP 09 (2026) 090, arXiv:2605.27641. [T34] | Robin propagator in BTZ; a starting point for finite temperature. | background | Thong refs |
