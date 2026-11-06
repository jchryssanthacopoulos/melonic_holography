# Melonic Holography: Robin Ladders and N=2 SYK

## Project Overview

This repository extends the construction of Thong, *Melonic ladder holography with a single free bulk field in AdS$_{d+1}$* (arXiv:2610.31989, `papers/thong_2026.pdf`) to $\mathcal{N}=2$ supersymmetric SYK.

**Thong's result.** Integrate out the centre of a CFT bilocal, and the pair Casimir acts on a relative cylinder $\mathbb{R}\times S^{d-1}$, with cylinder energy $\omega=h-d/2$. For a melonic CFT in the normalisable and stable range, the reduced ladder sum is the Robin boundary-to-boundary propagator of a single free scalar in rigid AdS$_{d+1}$. The constituent dimension fixes the mass and the ladder coupling fixes the Robin coefficient. This is a four-point statement, not an equivalence of theories.

**Our starting point.** `research/pdf/n2_susy_generalization_note.pdf` is an AI-generated draft proposing the extension: two scalars with a constant matrix Robin condition for the neutral channels, a Dirac response for the charged channels, and a fixed integral transform that reconstructs superconformal four-point functions. Treat it as a set of claims to test.

**Goals.**

1. Verify the draft against the primary literature (Peng–Spradlin–Volovich; Bulycheva).
2. Fill in the steps of Thong's derivation that the draft takes as given, notably a super-Casimir origin for the cylinder and the centre reduction of superblocks.
3. Understand what the auxiliary bulk system is: its validity range, field content, surviving symmetry, and the $H=1$ multiplet.

**Where things stand.** The living documents hold the current state; read these first:

- `docs/research_questions.md`: open questions and their status.
- `docs/derivations.md`: claim-by-claim status of the draft, notation, and project-derived results.
- `docs/todo.md`: the next tasks.
- `literature/synthesis.md`: the current understanding of the literature.

## Repository Structure

    papers/                  Primary literature (PDFs)
    literature/
      bibliography.md        All references, with reading status and metadata source
      synthesis.md           Current understanding across the literature
      notes/                 One note per reviewed paper
    docs/
      derivations.md         Derivations, notation, and status of claims
      research_questions.md  Open questions
      todo.md                Task list
    research/                Earlier, non-authoritative material
      pdf/                   Generated reports, including the N=2 draft
      notes/, tex/           Markdown notes and LaTeX sources (empty so far)
    scripts/                 Runnable experiments and checks
    src/                     Reusable research code (empty so far)
    results/
      data/                  Generated data and run logs
      figures/               Figures

The papers in `papers/` form the backbone of the literature review. They should be treated as primary sources for the theoretical framework and existing results.

## Conventions

Use Thong's notation wherever it applies, and the draft's additions otherwise. `docs/derivations.md` Sec. 1 has the dictionary, including collisions: $J$ is both a harmonic label and a coupling, Thong's $\mu(h)$ is unrelated to the draft's $\mu=\tfrac12-\delta$, and PSV's "s" and "a" labels swap between kernels and four-point functions. Cite equations as [T](3.11) for Thong arXiv v1, [N2](9) for the draft, and the tags in `literature/notes/` for other papers. The Maldacena–Stanford PDF is arXiv v1, so its numbers differ from those Thong cites.

Mark each substantial claim with its status: paper claim, draft claim, checked, derived here, or open. "Checked" means a formula passes a test. It does not mean the formula describes N=2 SYK.

## Research Philosophy

This is a physics research project, not merely a software project.

When working on the project:

- Prioritize mathematical and physical correctness over convenience.
- Clearly distinguish established results from conjectures.
- Do not silently assume results that have not been established.
- When deriving a result, show the derivation rather than simply stating it.
- Check limiting cases and dimensions whenever appropriate.
- Compare numerical results against analytical predictions whenever possible.
- Prefer reproducible numerical experiments.
- Record important computational choices.
- Do not modify existing results/data without explicitly explaining why.

When uncertain about a physics claim, investigate the relevant papers in `papers/` before making an assertion.

## Literature

The papers in `papers/` are the primary literature sources. For each important paper, maintain a corresponding note in:

    literature/notes/

Each paper note should contain:

- Full citation
- Main question
- Physical system
- Important definitions
- Main assumptions
- Main analytical results
- Important equations
- Numerical methods
- Relevant figures/results
- Limitations
- Relationship to our project

Whenever you review a paper:

1. Create or update its note in `literature/notes/`.
2. Update `literature/bibliography.md`, including its reading status.
3. Update `literature/synthesis.md`.
4. Update `docs/` where the paper changes a derivation, question, or task.

When citing a result from a paper, identify the paper and, when practical, the relevant equation/section/page. Do not attribute a result to a paper unless it is actually present there.

## Existing Research Material

The `research/` directory contains material produced during earlier stages of this project, including Claude-generated analyses, research notes, preliminary derivations, literature reviews, and other exploratory work.

This material is valuable context but is not automatically authoritative. In particular, Claude-generated documents may contain mistakes, unsupported claims, incorrect derivations, or interpretations that were subsequently revised. Treat them as prior research notes rather than as established facts.

When using material from `research/`:

- Preserve useful prior reasoning and context rather than unnecessarily recreating work.
- Distinguish clearly between claims supported by the primary literature, mathematical derivations performed within the project, numerical observations, and speculative hypotheses.
- When a claim from an existing research document is important to the current argument, verify it against the relevant primary paper, derivation, or numerical evidence.
- Do not propagate an unsupported claim simply because it appears in an existing project document.
- If an existing derivation appears questionable, identify and investigate the issue rather than silently treating it as correct.
- When the current understanding of the project changes, update the appropriate living documents rather than relying solely on historical notes.

The `research/` directory should therefore be viewed as the project's research history and working knowledge base, while `papers/` contains the primary literature and the current project documentation contains the best available statement of our present understanding.

In particular:

- `research/pdf/` contains previously generated research reports and analyses, including the N=2 draft. The draft cites a pre-arXiv version of Thong's paper, so its equation references do not match arXiv v1 (see `docs/derivations.md` Sec. 2.4).
- `research/notes/` contains earlier Markdown notes and exploratory material.
- `literature/synthesis.md` contains the current synthesis of the research.

Historical material should generally be preserved rather than deleted. If an old result is superseded, make that clear in the relevant document rather than erasing the historical record. Old historical results may also be saved
in the `results/data/` and `results/figures/` directories.

## Current Research Questions

See:

    docs/research_questions.md

Keep this document up to date as the research develops.

## Theoretical Work

Put substantial derivations and theoretical notes in:

    docs/derivations.md

When developing a derivation:

1. State the assumptions.
2. Define all variables.
3. Derive the result step by step.
4. Check special/limiting cases.
5. State clearly what has been proven versus conjectured.

## Code

Python is the primary programming language. Reusable research code belongs in `src/`, and runnable experiments and checks in `scripts/`.

Do not put substantial reusable code directly into notebooks. Code should be:

- reproducible
- reasonably modular
- documented
- numerically stable where possible

Use command-line arguments for experiment parameters rather than hard-coding them. For example:

    python scripts/run_experiment.py --q 3 --g 1.0 --parity 0 --dps 30

Currently the code is three check scripts. Run them after changing anything they cover.

- `scripts/check_robin_identities.py`: closed-form relations from Thong and the draft.
- `scripts/check_literature_inputs.py`: the draft against formulas transcribed from its source papers, with equation numbers. It imports from the first script.
- `scripts/check_supercasimir.py`: the $su(1,1|1)$ generators and pair Casimir, as operators on a Grassmann algebra.

The project has no `.venv` yet: these use system `python3` with `mpmath` and `sympy`, and SciPy in the current Anaconda install fails to import.

Code produced with the assistance of Claude should be treated as a starting point for further development, not as automatically correct or authoritative.

Before writing new code:

- Inspect the existing codebase and determine whether relevant functionality already exists.
- Prefer extending, refactoring, or validating existing code over unnecessarily duplicating it.
- Understand the purpose and assumptions of existing code before modifying it.
- Do not silently discard or replace existing implementations.
- Treat existing Claude-generated code as potentially useful but unverified. Check mathematical formulas, physical assumptions, numerical methods, normalization conventions, and edge cases before relying on important results.
- If existing code appears incorrect, identify the issue explicitly and explain the proposed correction.
- Preserve working functionality unless there is a clear reason to change it.
When making substantial changes, explain what was changed and why.

## Numerical Experiments

Each significant numerical experiment should record:

- parameters
- random seed(s), when applicable
- number of disorder realizations
- system size
- theoretical prediction being tested
- output files
- relevant code version/commit

### Memory safety (mandatory)

**Every intensive process must run under the memory watchdog.** This machine has 32 GB of RAM; a campaign on
2026-09-28 was OOM-killed when two large jobs ran concurrently. The project budget is **30 GB total** across all
our running jobs.

**`scripts/memwatch.py` is not in this repository yet.** Add it (see `docs/todo.md`) before launching any intensive job. The description below is how it is expected to behave.

- Before launching, estimate the peak memory, e.g. with a `--dry-run` option in the experiment script. Do not launch a job whose estimate does not fit.
- Launch through `scripts/memwatch.py`, with a per-job cap somewhat above the estimate:

      .venv/bin/python scripts/memwatch.py --limit-gb 12 --log results/data/<run>/memwatch.jsonl -- <command>

  memwatch measures the whole process tree (macOS physical footprint, which counts compressed memory), kills the
  job if it exceeds its cap, if all watched jobs together exceed the 30 GB budget, or if the system runs low. Each
  job reserves its cap: a launch that would push the reservations of running jobs over 30 GB is refused (or
  queued with `--wait`). So concurrent jobs' caps must sum to ≤ 30 GB. `memwatch.py --status` lists the jobs.
- "Intensive" means anything that may exceed ~2 GB or run longer than a few minutes, including parallel groups
  of small jobs. When in doubt, use the watchdog.
- Prefer algorithms whose memory is bounded (bounded caches, partial eigendecompositions, streaming) over ones
  that hold several dense D×D matrices at once. Record the measured peak (memwatch prints it) with the results.

Results should go into:

    results/

Figures should go into:

    results/figures/

Raw/generated data should go into:

    results/data/

Do not overwrite previous experimental results without a good reason.

## Figures

Figures should be publication-quality whenever practical. Every figure should have:

- labeled axes
- units where appropriate
- legends when needed
- sensible scaling
- readable fonts
- a clear caption or accompanying description

Whenever comparing numerical data to an analytical prediction, plot the analytical prediction explicitly.

## Git

This directory is not a git repository yet (see `docs/todo.md`). Until it is, record the date in place of a commit in experiment records.

Make commits at meaningful milestones. Do not commit:

- huge generated datasets
- temporary files
- Python cache files
- credentials
- unnecessary notebook outputs

Before making a substantial change, inspect the existing code and understand the current implementation rather than rewriting things unnecessarily.

## Claude's Role

Act as a research collaborator and computational physicist.

You should:

- read relevant papers before proposing literature-based claims
- help formulate hypotheses
- derive equations
- design numerical experiments
- implement simulations
- analyze results
- identify inconsistencies
- suggest useful tests
- maintain reproducibility

Do not merely agree with my hypotheses. If something appears incorrect, say so explicitly and explain why. When there are multiple plausible interpretations, state them explicitly and determine which is consistent with the papers and existing project.

## Important Rule

Never fabricate a citation, equation, numerical result, or claim about a paper. If the necessary information is not available, say so and investigate the available project materials before proceeding.
