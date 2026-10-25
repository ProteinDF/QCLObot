#!/usr/bin/env bash
# scripts/ci-local.sh: Run CI checks locally
set -uo pipefail

# Always run from the repository root (CI_LOCAL_ROOT checks another worktree with
# this copy of the script; devtool/review.sh uses it)
REPO_ROOT="$(cd "${CI_LOCAL_ROOT:-$(dirname "${BASH_SOURCE[0]}")/..}" && pwd)"

# Prevent Python bytecode caching (__pycache__) in repository
export PYTHONDONTWRITEBYTECODE=1

PYTHON="${PYTHON:-python3}"

# Jobs run on a copy of the working tree in a temporary directory outside the
# repository, so test outputs and Sphinx builds never land in the worktree.
TMP_ROOT="$(mktemp -d -t qclobot-ci-XXXXXX)"
cleanup() {
    rm -rf "$TMP_ROOT"
}
trap cleanup EXIT

ALL_DEFAULT_JOBS=("unit" "docs")
REQUESTED_JOBS=()

show_help() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] [JOB...]

Run QCLObot CI checks locally.
If no jobs are specified, all default jobs are executed:
  ${ALL_DEFAULT_JOBS[*]}

Options:
  -h, --help       Show this help message and exit

Available Jobs:
  unit             Run unittest (tests/test_*.py) on a copy of the working tree.
                   Requires AMBERHOME (with bin/tleap, bin/sander, bin/ambpdb,
                   bin/reduce) and an importable proteindf_bridge.
  docs             Build English and Japanese Sphinx documentation in an isolated venv
  all              Run all default jobs (${ALL_DEFAULT_JOBS[*]})

Environment Variables:
  PYTHON           Python interpreter to use (default: python3)
  AMBERHOME        AmberTools installation (required for the unit job)
  CI_LOCAL_ROOT    Repository root to check (default: the one containing this script)

Prerequisites:
  - Python 3.9+ with venv support
  - AmberTools (tleap, sander, ambpdb, reduce) for the unit job
  - ProteinDF_bridge installed in \$PYTHON (imported by the tests)

Failure Handling:
  All requested jobs are executed even if an earlier job fails. A summary of
  SUCCESS/FAILED and duration per job is printed at the end. If any job failed,
  the script exits with status code 1. A missing prerequisite is a failure,
  not a skip.
EOF
}

# Copy the working tree (tracked and untracked, minus ignored files) into a
# temporary directory, so that the checks run outside the repository.
copy_tree() {
    local dest="$TMP_ROOT/src"
    mkdir -p "$dest"
    (cd "$REPO_ROOT" && git ls-files -co --exclude-standard -z) \
        | (cd "$REPO_ROOT" && tar --null -T - -cf -) \
        | tar -xf - -C "$dest"
    echo "$dest"
}

# ---------------------------------------------------------------------------
# Jobs
# ---------------------------------------------------------------------------

job_unit() {
    echo "[unit] preflight"
    if [ -z "${AMBERHOME:-}" ]; then
        echo "[unit] error: AMBERHOME is not set" >&2
        return 1
    fi
    local tool
    for tool in tleap sander ambpdb reduce; do
        if [ ! -x "$AMBERHOME/bin/$tool" ]; then
            echo "[unit] error: missing executable: $AMBERHOME/bin/$tool" >&2
            return 1
        fi
    done
    if ! "$PYTHON" -c "import proteindf_bridge" 2>/dev/null; then
        echo "[unit] error: proteindf_bridge is not importable with $PYTHON" >&2
        return 1
    fi

    local src
    src="$(copy_tree)" || return 1
    cd "$src" || return 1
    echo "[unit] python -m unittest discover -s tests -t ."
    "$PYTHON" -m unittest discover -s tests -t . -p 'test_*.py' 2>&1
}

job_docs() {
    local venv="$TMP_ROOT/venv" src
    src="$(copy_tree)" || return 1

    echo "[docs] creating isolated venv"
    "$PYTHON" -m venv "$venv" || return 1
    "$venv/bin/python" -m pip install -q --upgrade pip || return 1
    "$venv/bin/python" -m pip install -q -r "$src/doc/requirements.txt" || return 1

    cd "$src/doc" || return 1
    echo "[docs] English HTML"
    make html SPHINXBUILD="$venv/bin/sphinx-build" || return 1
    echo "[docs] Japanese HTML"
    make -e SPHINXOPTS="-D language='ja'" BUILDDIR=_build_ja \
        SPHINXBUILD="$venv/bin/sphinx-build" html || return 1
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help) show_help; exit 0 ;;
        all) REQUESTED_JOBS+=("${ALL_DEFAULT_JOBS[@]}") ;;
        unit|docs) REQUESTED_JOBS+=("$1") ;;
        *) echo "error: unknown option or job: $1" >&2; show_help >&2; exit 2 ;;
    esac
    shift
done
if [ ${#REQUESTED_JOBS[@]} -eq 0 ]; then
    REQUESTED_JOBS=("${ALL_DEFAULT_JOBS[@]}")
fi

declare -A RESULT DURATION
overall=0
for job in "${REQUESTED_JOBS[@]}"; do
    start=$(date +%s)
    echo "=============== job: $job ==============="
    (job_"$job")
    status=$?
    DURATION[$job]=$(( $(date +%s) - start ))
    if [ $status -eq 0 ]; then
        RESULT[$job]="SUCCESS"
    else
        RESULT[$job]="FAILED"
        overall=1
    fi
    echo
done

echo "=============== CI Local Summary ==============="
printf '%-10s %-10s %s\n' JOB RESULT DURATION
for job in "${REQUESTED_JOBS[@]}"; do
    printf '%-10s %-10s %ss\n' "$job" "${RESULT[$job]}" "${DURATION[$job]}"
done
if [ $overall -eq 0 ]; then
    echo "Result: all requested jobs passed"
else
    echo "Result: one or more jobs FAILED"
fi
exit $overall
