# Melonic Holography: Robin Ladders and N=2 SYK

A research project extending the construction of [Thong (2026)](https://arxiv.org/abs/2610.31989) to $\mathcal{N}=2$ supersymmetric SYK.

Thong shows that, once the pair centre is integrated out, the ladder sum of a melonic CFT$_d$ is the Robin boundary-to-boundary propagator of a single free scalar in rigid AdS$_{d+1}$. The ladder coupling sits entirely in the boundary condition. That result covers one ladder with one constituent dimension. $\mathcal{N}=2$ SYK has mixed fermion–boson channels, charged channels and superconformal multiplets, so it needs more.

An AI-generated draft, [`n2_susy_generalization_note.pdf`](research/pdf/n2_susy_generalization_note.pdf), proposes two scalars with a matrix Robin condition, plus a Dirac response and a superconformal reconstruction. The goal here is to verify that draft against the primary literature, fill the gaps in its derivation, and work out what the auxiliary bulk system means.

## Status

Early, as of November 2026. Ten further papers have been reviewed, some in part. The main results so far:

- **The draft's physical input agrees with its sources,** and the auxiliary covariance reproduces the literature's centre-integrated four-point function at finite coupling, where tested.
- **The relative cylinder follows from the super-Casimir,** with a purely bosonic relative space and $\omega=H$.
- **The apparent contradiction at the $H=1$ multiplet is resolved.**
- **The super-centre integral of the four-point function is exactly the Dirac propagator,** with one exponential per superblock.
- **The bulk boundary condition breaks $su(1,1|1)$ to $u(1|1)$:** besides the AdS isometries $E_\pm$, two of the four Killing supercharges are lost, in each embedding of the algebra on the fields.

The write-up is [`research/pdf/n2_dirac_ladders.pdf`](research/pdf/n2_dirac_ladders.pdf), and its [audit](research/notes/n2_dirac_ladders_audit.md) says how far to trust each claim. The remaining gaps are listed in the paper's App. B and in [the derivations](docs/derivations.md). These project results have not been externally reviewed or searched for in prior work. See [the derivations](docs/derivations.md), [the open questions](docs/research_questions.md) and [the task list](docs/todo.md).

## Layout

| Path | Contents |
|---|---|
| [`papers/`](papers/) | Primary literature (PDFs) |
| [`literature/`](literature/) | [Bibliography](literature/bibliography.md), [synthesis](literature/synthesis.md), and per-paper [notes](literature/notes/) |
| [`docs/`](docs/) | [Derivations and claim status](docs/derivations.md), [research questions](docs/research_questions.md), [to-do list](docs/todo.md) |
| [`research/`](research/) | AI-generated drafts and reports, including the working paper ([`tex/`](research/tex/), [`pdf/`](research/pdf/)) and its audit ([`notes/`](research/notes/)). Useful, but not authoritative. |
| [`scripts/`](scripts/) | Runnable experiments and checks |
| [`src/`](src/) | Reusable research code (empty so far) |
| [`results/`](results/) | Generated [data](results/data/) and [figures](results/figures/) |

## Running the checks

The checks need `mpmath` and, for the last, `sympy`:

```bash
python3 scripts/check_robin_identities.py --dps 30   # Thong and the draft: internal identities
python3 scripts/check_literature_inputs.py --dps 25  # the draft against its source papers
python3 scripts/check_supercasimir.py                # the su(1,1|1) pair Casimir
python3 scripts/check_dirac_reduction.py --dps 20    # the reduced kernel as a Dirac response
python3 scripts/check_bulk_susy_algebra.py           # Killing spinors and closure (about 2 minutes)
python3 scripts/check_bulk_susy_modes.py --dps 20    # supercharges on normal modes and boundary conditions
python3 scripts/check_ward_and_vacuum.py --dps 20    # Ward identity, threshold, figure data (with --write)
python3 scripts/check_reconstruction.py --dps 18     # the four-point reconstruction at finite g
python3 scripts/check_super_centre.py --dps 20       # the super-centre integral (about 1 minute)
python3 scripts/audit_n2_paper.py --dps 20           # independent audit of the paper (about 2 minutes)
```

Each exits non-zero if a check fails. Recorded runs are in [`results/data/`](results/data/).

## Working with Claude

This project uses Claude Code. See [CLAUDE.md](CLAUDE.md) for its instructions.
