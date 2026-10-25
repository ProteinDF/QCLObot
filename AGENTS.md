# Rules for implementation agents (agy / GitHub Copilot)

> 日本語版: [AGENTS.ja.md](AGENTS.ja.md)

In this repository, implementation is done by agy or GitHub Copilot, review by Claude, and merge approval by the user.
The overall flow is in `CONTRIBUTING.md`. This file lists only what implementation agents must follow.

## Working location and branch (MUST)

- Work in the **dedicated worktree** of the branch given in the task file (the `TASK_*.md` path that `devtool/delegate.sh` passes). The worktree is created by `devtool/flow.sh start <branch>` (or by `devtool/delegate.sh`), next to the `develop` worktree, with `/` in the branch name replaced by `-`.
- **Do not edit files in the `develop` worktree or in other worktrees.** Read the task file but do not change it (Claude appends review results).
- The task file is private. Do not copy its text (for example the analysis of a bug) into commit messages or code comments.
- Do not write the task file name (`TASK_<name>`) in commit messages, code comments, docstrings, or documents. In test descriptions, say what is checked.
- Do not write local directory paths (for example a working location under the home directory) in files that are committed or in commit messages. Read local locations from environment variables or from a git-ignored configuration file.
- Do not commit to any other branch. Do not merge into `develop` or `main`, do not push, do not create tags or releases, do not delete branches or worktrees, and do not use `git stash`.
- Check the location and branch with `git status` and `git branch --show-current` before you start work.
- Put temporary files and temporary virtual environments outside the repository (for example with `mktemp -d`). Do not create `.venv` or `dist/` in the repository.

## Implementation

- Follow the task file's "targets" and "definition of done". If a design decision is needed that the task file does not cover, do not decide it yourself: write it in the "decisions needed" part of the completion report and stop.
- Do not make out-of-scope changes (unrelated refactoring, formatting-only diffs, mass file moves).
- Match the style of the existing code. Change only the lines you need to change; do not reformat whole existing files.
- Do not swallow errors (for example by silently replacing them with default values). Do not hide warnings.
- Do not weaken tests: do not remove assertions silently, and do not add no-op methods or duplicate aliases to the public API.
- Commits should be logical units. Follow the existing history for messages: `feat: ...`, `fix: ...`, `test: ...`, `docs: ...`, `chore: ...` (in English).
- Documentation changes must include the Japanese translation in the same commit (see `CONTRIBUTING.md`).

## Context-efficient reading and running

- Do not read large files in full. Locate symbols with `rg -n '<symbol>' <dir>` and read only the needed range with `sed -n '<start>,<end>p' <file>`.
- Save build and test output to a file and look only at the tail or the failing parts.
- Do not re-read the same file repeatedly. Read the task file once, at the start.
- Do not print progress messages while waiting.

## Verification before completion

Run this in the worktree and confirm it succeeds (see `CONTRIBUTING.md`, "Running tests and CI"):

```bash
scripts/ci-local.sh
```

The `unit` job needs AmberTools (`AMBERHOME`) and `proteindf_bridge`. If they are not available, `unit` fails by design. Say so in the report and do not treat that as a pass.

## Completion report

Finish with a report in this form. **Do not copy numbers from memory: paste the matching lines of the actual output.**
Keep it **short**: paste only the summary and result lines (the `CI Local Summary` section and test counts), not full build or test logs. The reviewer reruns the checks.
Confirm test names, file names and line numbers by opening the actual files before you write them. Do not state unverified guesses as facts.
For each item in the definition of done, write one or two lines on what you did.

```markdown
## Completion report: <branch>

### Commits
<output of git log --oneline develop..HEAD>

### What was done
<for each item of the definition of done, what was done; match item numbers>

### Verification
<the CI Local Summary section of scripts/ci-local.sh, pasted as is>
<result lines for any extra checks the task file asks for>

### Decisions needed / open items
<"None" if there are none>
```
