# Melonic Holography: Robin Ladders and N=2 SYK

A research project extending the construction of [Thong (2026)](https://arxiv.org/abs/2610.31989) to $\mathcal{N}=2$ supersymmetric SYK.

Thong shows that, once the pair centre is integrated out, the ladder sum of a melonic CFT$_d$ is the Robin boundary-to-boundary propagator of a single free scalar in rigid AdS$_{d+1}$. The ladder coupling sits entirely in the boundary condition. That result covers one ladder with one constituent dimension. $\mathcal{N}=2$ SYK has mixed fermion–boson channels, charged channels and superconformal multiplets, so it needs more.

An AI-generated draft, [`n2_susy_generalization_note.pdf`](research/pdf/n2_susy_generalization_note.pdf), proposes two scalars with a matrix Robin condition, plus a Dirac response and a superconformal reconstruction. The goal here is to verify that draft against the primary literature, fill the gaps in its derivation, and work out what the auxiliary bulk system means.

## Status

Early, as of November 2026. Ten further papers have been reviewed, some in part. The main results so far:

- **The draft's physical input agrees with its sources,** and the auxiliary covariance reproduces the literature's centre-integrated four-point function at finite coupling, where tested.
- **The relative cylinder follows from the super-Casimir,** with a purely bosonic relative space and $\omega=H$.
- **The apparent contradiction at the $H=1$ multiplet is resolved.**

The main remaining gap is the super-centre integral of a superblock. These project results have not been externally reviewed or searched for in prior work. See [the derivations](docs/derivations.md), [the open questions](docs/research_questions.md) and [the task list](docs/todo.md).

## Layout

| Path | Contents |
|---|---|
| [`papers/`](papers/) | Primary literature (PDFs) |
| [`literature/`](literature/) | [Bibliography](literature/bibliography.md), [synthesis](literature/synthesis.md), and per-paper [notes](literature/notes/) |
| [`docs/`](docs/) | [Derivations and claim status](docs/derivations.md), [research questions](docs/research_questions.md), [to-do list](docs/todo.md) |
| [`research/`](research/) | Earlier material, including AI-generated drafts. Useful, but not authoritative. |
| [`scripts/`](scripts/) | Runnable experiments and checks |
| [`src/`](src/) | Reusable research code (empty so far) |
| [`results/`](results/) | Generated [data](results/data/) and [figures](results/figures/) |

## Running the checks

The checks need `mpmath` and, for the last, `sympy`:

```bash
python3 scripts/check_robin_identities.py --dps 30   # Thong and the draft: internal identities
python3 scripts/check_literature_inputs.py --dps 25  # the draft against its source papers
python3 scripts/check_supercasimir.py                # the su(1,1|1) pair Casimir
```

Each exits non-zero if a check fails. Recorded runs are in [`results/data/`](results/data/).

## Working with Claude

This project uses Claude Code. See [CLAUDE.md](CLAUDE.md) for its instructions.
