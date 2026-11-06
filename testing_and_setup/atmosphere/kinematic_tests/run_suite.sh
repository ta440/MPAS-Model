#!/bin/bash
#
# Run a suite of kinematic transport tests. See README.md.
#
# Usage: run_suite.sh <quick|convergence|scaling> -d <meshes_dir> [-b <exe_dir>] [-n <ranks>] [-k <rk_order>]
#                    [-L <on|off>]
#   <meshes_dir> holds one subdirectory per mesh (x1.10242/, x1.40962/, x1.163842/).
#
HERE=$(cd "$(dirname "$0")" && pwd)
SUITE=$1; shift || true
MESHES=; EXTRA=()
while getopts "d:b:n:k:L:" opt; do
   case $opt in
      d) MESHES=$OPTARG ;;
      b) EXTRA+=(-b "$OPTARG") ;;
      n) EXTRA+=(-n "$OPTARG") ;;
      k) EXTRA+=(-k "$OPTARG") ;;
      L) EXTRA+=(-L "$OPTARG") ;;
      *) sed -n '2,7p' "$0"; exit 2 ;;
   esac
done
[ -n "$SUITE" ] && [ -n "$MESHES" ] || { sed -n '2,7p' "$0"; exit 2; }

FAILED=()
run() {   # run <mesh> <extra options...>
   local mesh=$1; shift
   "$HERE/run_kinematic_test.sh" -c lauritzen_nondiv -m $mesh -d "$MESHES/$mesh" "${EXTRA[@]}" "$@" \
      || FAILED+=("$mesh $*")
}

case $SUITE in
   quick)
      run x1.10242 ;;
   convergence)
      for mesh in x1.10242 x1.40962 x1.163842; do run $mesh; done ;;
   scaling)
      # About 11 tracers in a default MPAS run; up to 200 for chemistry
      for n in 4 11 50 100 200; do
         run x1.40962 -N $n -r 1_00:00:00 -i 1_00:00:00 -o "$PWD/runs/scaling_x1.40962_N${n}"
      done ;;
   *)
      echo "Unknown suite '$SUITE' (quick, convergence or scaling)"; exit 2 ;;
esac

if [ ${#FAILED[@]} -eq 0 ]; then
   echo "Suite '$SUITE': all runs passed"
else
   echo "Suite '$SUITE': ${#FAILED[@]} run(s) failed:"; printf '   %s\n' "${FAILED[@]}"; exit 1
fi
