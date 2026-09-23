#!/usr/bin/env bash
# Automated pre-push checks: pre-commit hooks, nf-test suite, and nf-core
# lint. Fail-fast (set -e) so the first broken step stops the script.
#
# Note on dependencies: the dry-run and lint steps need nothing but the
# tools themselves, but the real (non-dry-run) nf-test suite still spins up
# Docker containers, and its pipeline-level tests (tests/default.nf.test)
# need the S3 test dataset synced locally (data-test/ -- see CLAUDE.md's
# "Test dataset" section). This script does NOT include the ad hoc, manually
# -inspected `nextflow run` smoke tests -- see scripts/run-smoke-tests.sh for
# those.
#
# Usage: scripts/run-test-suite.sh

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"
export NXF_FILE_ROOT="$PWD"

require() {
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found on PATH." >&2; exit 1; }
}
require pre-commit
require nf-test
require nf-core

step() { echo; echo "==> $1"; }

step "[1/4] pre-commit (prettier, editorconfig-checker, ...)"
pre-commit run --all-files

# Remember that running 'nf-test test' essentially discovers and run .nf.test files
step "[2/4] nf-test dry-run (syntax check only, no execution)"
nf-test test --dryRun

step "[3/4] nf-test full suite (--profile test,docker)"
nf-test test --profile test,docker

step "[4/4] nf-core pipelines lint --release"
nf-core pipelines lint --release

echo
echo "All checks passed."
