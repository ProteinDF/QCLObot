#!/usr/bin/env bash
# Reviewer-side verification of a working branch (run by Claude, from any worktree).
# Runs scripts/ci-local.sh against the branch's worktree and saves the full output
# for the review record.
#
# Usage: devtool/review.sh [--verbose] <branch> [ci-local.sh job...]
#   jobs      passed to ci-local.sh (default: all default jobs)
#   --verbose print the whole output (default: header, then the summary only)
#
# The reviewer's (develop) copy of ci-local.sh is used, with CI_LOCAL_ROOT set to the
# branch's worktree, so an edited copy on the branch cannot weaken the checks.
# The full output is always saved under <git-common-dir>/review-logs/.
set -uo pipefail

here=$(cd "$(dirname "$0")" && pwd)
verbose=0
while (($#)); do
    case "$1" in
        --verbose) verbose=1; shift ;;
        -*) sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
        *) break ;;
    esac
done
(($# >= 1)) || { sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
branch=$1; shift
jobs=("$@")

wt=$("$here/flow.sh" path "$branch") || exit 1
logdir="$(git rev-parse --path-format=absolute --git-common-dir)/review-logs"
mkdir -p "$logdir"
log="$logdir/$(date +%Y%m%d-%H%M%S)_${branch//\//-}.log"

run() {
    echo "branch:   $branch"
    echo "worktree: $wt"
    echo
    echo "--- commits (develop..$branch)"
    git -C "$wt" log --oneline "develop..$branch"
    echo
    echo "--- diffstat"
    git -C "$wt" diff --stat "develop...$branch"
    echo
    echo "--- worktree status (uncommitted / untracked)"
    git -C "$wt" status --short
    echo

    echo "--- ci-local.sh ${jobs[*]:-(default jobs)}"
    CI_LOCAL_ROOT="$wt" "$here/../scripts/ci-local.sh" "${jobs[@]}"
    ci_status=$?

    echo
    echo "--- worktree status after ci-local.sh (must not grow)"
    git -C "$wt" status --short
    echo
    echo "=== review: ci-local.sh exit $ci_status"
    ((ci_status == 0))
}

if ((verbose)); then
    run 2>&1 | tee "$log"
    status=${PIPESTATUS[0]}
else
    run >"$log" 2>&1
    status=$?
    # header (up to the diffstat), the ci-local.sh summary and the verdict
    sed -n '1,/^--- worktree status/p' "$log" | sed '$d'
    sed -n '/^CI Local Summary/,/^Result:/p' "$log"
    sed -n '/^--- worktree status (/,/^--- ci-local.sh/p' "$log" | sed '1d;$d' | grep . | sed 's/^/uncommitted: /'
    sed -n '/^--- worktree status after/,/^=== review/p' "$log" | sed '1d;$d' | grep . | sed 's/^/after ci: /'
    grep '^=== review' "$log"
fi
echo "log: $log"
exit "$status"
