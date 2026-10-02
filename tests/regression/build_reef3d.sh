#!/usr/bin/env bash
# Build a REEF3D binary for the regression suite.
#
#   build_reef3d.sh -o <outdir> [-r <git rev>] [-j <jobs>] [-H <hypre dir>] [-I <hypre include dir>]
#
#   -o  output directory; the binary is <outdir>/REEF3D, objects in <outdir>/obj
#       (incremental: re-running only recompiles changed files)
#   -r  build a git revision (branch, tag, commit) instead of the working tree;
#       it is checked out as a detached git worktree in <outdir>/src (the working tree is not touched;
#       remove later with 'git worktree remove <outdir>/src')
#   -j  parallel jobs (default: number of cores)
#   -H  hypre prefix (default: $HYPRE_DIR or /usr/local/hypre)
#   -I  hypre include dir (default: <hypre prefix>/include)
#
# Flags are fixed so that two binaries built from different sources can be compared bitwise:
#   -O2 -ffp-contract=off   (no FMA contraction, no -march=native, no -ffast-math, no LTO)
# Use the same compiler and MPI for both binaries of an A/B comparison.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUT=""
REV=""
JOBS="$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
HYPRE="${HYPRE_DIR:-/usr/local/hypre}"
HYPRE_INC=""

while getopts "o:r:j:H:I:h" opt; do
    case $opt in
        o) OUT="$OPTARG" ;;
        r) REV="$OPTARG" ;;
        j) JOBS="$OPTARG" ;;
        H) HYPRE="$OPTARG" ;;
        I) HYPRE_INC="$OPTARG" ;;
        *) sed -n '2,20p' "$0"; exit 1 ;;
    esac
done
[ -z "$OUT" ] && { sed -n '2,20p' "$0"; exit 1; }
HYPRE_INC="${HYPRE_INC:-$HYPRE/include}"

mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"

if [ -n "$REV" ]; then
    # detached git worktree of the revision: the main working tree is not touched, and switching
    # to another revision later only touches changed files (incremental rebuild)
    SRC="$OUT/src"
    COMMIT="$(git -C "$REPO" rev-parse --short=7 "$REV")"
    BRANCH="$REV"
    if [ ! -e "$SRC/.git" ]; then
        git -C "$REPO" worktree add --detach "$SRC" "$REV" >/dev/null
    else
        git -C "$SRC" checkout -q --detach "$REV"
    fi
else
    SRC="$REPO"
    COMMIT="$(git -C "$REPO" rev-parse --short=7 HEAD)$(git -C "$REPO" diff --quiet HEAD -- src || echo -dirty)"
    BRANCH="$(git -C "$REPO" rev-parse --abbrev-ref HEAD)"
fi

if [ ! -f "$SRC/src/regression_dump.cpp" ]; then
    echo "warning: $SRC has no src/regression_dump.cpp - the binary will not write regression dumps" >&2
fi

echo "building REEF3D $BRANCH@$COMMIT from $SRC -> $OUT/REEF3D (jobs=$JOBS)"
make -C "$SRC" -j "$JOBS" all \
    OBJ_DIR="$OUT/obj" APP_DIR="$OUT" HYPRE_DIR="$HYPRE" \
    GIT_BRANCH="$BRANCH" GIT_COMMIT="$COMMIT" GIT_DIRTY="" GIT_VERSION="$COMMIT" \
    CXXFLAGS="-std=c++20 -O2 -ffp-contract=off -w -DVERSION=\\\"$COMMIT\\\" -DBRANCH=\\\"$BRANCH\\\" -DBUILD=\\\"regression\\\"" \
    INCLUDE="-I $HYPRE_INC -I ThirdParty/eigen-5.0.0 -DEIGEN_MPL2_ONLY" \
    LDFLAGS="-L $HYPRE/lib -lHYPRE"

echo "$BRANCH@$COMMIT" > "$OUT/REEF3D.version"
echo "done: $OUT/REEF3D"
