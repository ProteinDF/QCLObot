# Contributing

> 日本語版: [CONTRIBUTING.ja.md](CONTRIBUTING.ja.md)

This file is the English original. When you change it, update `CONTRIBUTING.ja.md` in the same commit.

## Branching model (GitFlow)

This project follows a GitFlow-style branching model.

- **`main`**: the latest released state. Nothing is committed here directly. It only receives merges from `release/*` and `hotfix/*` branches, each merge tagged with the version number (for example `2026.10.0`, without a `v` prefix).
- **`develop`**: the integration branch for ongoing development. Day-to-day work targets this branch.
- **`feature/*`, `fix/*`, `chore/*`**: branch off `develop` for a unit of work and merge back into `develop` after review.
- **`release/X.Y.Z`**: branched from `develop` when preparing a release. Final adjustments (version bump, etc.) are made here, then it merges into both `main` (tagged) and `develop`.
- **`hotfix/*`**: branched from `main` for urgent fixes to a released version. Merges into both `main` (tagged) and `develop`.

Each working branch has its own git worktree, created next to the `develop` worktree (`devtool/flow.sh start <branch>`).

## Versioning

Versions follow `YYYY.M.PATCH` (calendar versioning). The first release in a given year and month is `PATCH=0` (for example `2026.10.0`); later releases in the same month increment `PATCH`.

The version is stored in `qclobot/__version__.py` as the `__version_info__` tuple, which `setup.cfg` reads. `devtool/flow.sh start release/<version>` and `hotfix/<version>` update it.

## Running tests and CI

Continuous integration runs locally with `scripts/ci-local.sh`.

### Prerequisites

- Python 3.9 or later, with the `venv` module (set `PYTHON` to choose the interpreter; the default is `python3`).
- ProteinDF_bridge installed in that interpreter. The tests import `proteindf_bridge`.
- AmberTools installed, with `AMBERHOME` pointing to it. The `unit` job needs `$AMBERHOME/bin/tleap`, `sander`, `ambpdb` and `reduce`, which the tests call through `AMBERHOME`.

### Running local CI

```bash
scripts/ci-local.sh              # all default jobs: unit, docs
scripts/ci-local.sh docs         # one job
scripts/ci-local.sh --help       # usage
```

Jobs:

- `unit`: runs `unittest` on `tests/test_*.py` against a copy of the working tree in a temporary directory, so test outputs do not appear in the worktree. It fails first if `AMBERHOME` or one of the AmberTools executables is missing, or if `proteindf_bridge` cannot be imported. A missing prerequisite is a failure, not a skip.
- `docs`: builds the English and the Japanese Sphinx HTML in an isolated virtual environment, from a copy of the working tree.

## Documentation (English and Japanese)

English is the original. Every documentation change must include its Japanese translation in the same commit:

- Sphinx sources (`doc/source/*.rst`): after editing English sources, update the Japanese catalogs under `doc/source/locale/ja/LC_MESSAGES/` (regenerate the `.pot` files and run `sphinx-intl update`, then translate the new messages).
- Markdown documents (`*.md`): keep the Japanese file `*.ja.md` next to the English file and update both.

The `docs` job builds both languages and fails if either build fails.

## API reference (`doc/source`)

After adding, removing or renaming a module under `qclobot/`, update the Sphinx API pages in `doc/source/` and rebuild the documentation with the `docs` job.

## Development with AI agents

| Who | Role | Rules |
|---|---|---|
| User | direction and design decisions, approving merges, pushing | |
| Claude | writing task files, reviewing, merging after approval | `CLAUDE.md` |
| agy / GitHub Copilot | implementation, in the working branch's own worktree | `AGENTS.md` |

1. **Task file**: Claude writes `TASK_<name>.md` with the working branch (`**Branch**:` line), the targets and the definition of done. Task files are private. They live in the `ProteinDF_private` repository under `qclobot/tasks/`, linked from the `develop` worktree as `private/` (ignored by git).
2. **Implementation**: hand the task to an agent. The branch and its worktree are created next to the `develop` worktree if they are missing.
   ```bash
   devtool/delegate.sh private/tasks/TASK_<name>.md            # agy (non-interactive, auto-approve)
   devtool/delegate.sh -a copilot private/tasks/TASK_<name>.md # GitHub Copilot CLI
   devtool/delegate.sh -i private/tasks/TASK_<name>.md         # interactive
   ```
   Output is saved under `.git/agent-logs/`. Extra directories an agent may access (for example a local ProteinDF_bridge checkout or the AmberTools installation) are listed one per line in `$(git rev-parse --git-common-dir)/delegate.conf`. When an agent stops on its usage limit, the script exits with 75; resume with `-c`, or wait with `-W`.
3. **Review**: Claude runs `devtool/review.sh <branch>`, which runs the `develop` copy of `scripts/ci-local.sh` against the branch's worktree, and records the result in the task file.
4. **Merge**: after the user approves, `devtool/flow.sh finish <branch>`. Releases use `devtool/flow.sh start release/<YYYY.M.PATCH>` (bumps the version) and `devtool/flow.sh finish release/<YYYY.M.PATCH>` (merges into `main`, tags, merges back into `develop`).

`devtool/tasks.sh` lists unmerged tasks and working branches (`--all` includes merged tasks). Nothing is pushed by these scripts.

## Public and private information

Public files (tracked files, commit messages, and documentation) must not contain local directory paths, task file names, or the contents of the private repository. Keep local paths in environment variables or in `$(git rev-parse --git-common-dir)/delegate.conf`.
