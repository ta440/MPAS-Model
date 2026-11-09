#!/bin/bash
#
# Run a suite of kinematic transport tests. See README.md.
#
# Usage: run_suite.sh <quick|convergence|scaling|ranks> -d <meshes_dir> [-b <exe_dir>] [-n <ranks>] [-k <rk_order>]
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
run() {   # run <case> <mesh> <extra options...>
   local case=$1 mesh=$2; shift 2
   "$HERE/run_kinematic_test.sh" -c $case -m $mesh -d "$MESHES/$mesh" "${EXTRA[@]}" "$@" \
      || FAILED+=("$case $mesh $*")
}

case $SUITE in
   quick)
      run lauritzen_nondiv x1.10242
      run lauritzen_div x1.10242 ;;
   convergence)
      for case in lauritzen_nondiv lauritzen_div; do
         for mesh in x1.10242 x1.40962 x1.163842; do run $case $mesh; done
      done ;;
   scaling)
      # Whole blocks of the five tracer shapes; N = 10 stands in for a default MPAS run
      # (about 11 tracers), N = 200 is the upper bound for chemistry
      for n in 5 10 20 50 100 200; do
         run lauritzen_nondiv x1.40962 -N $n -r 1_00:00:00 -i 1_00:00:00 -o "$PWD/runs/scaling_x1.40962_N${n}"
      done ;;
   ranks)
      # Results must not depend on the decomposition: 1 and 4 ranks, all fields bitwise equal
      for case in lauritzen_nondiv lauritzen_div; do
         for n in 1 4; do
            run $case x1.10242 -n $n -o "$PWD/runs/ranks_${case}_np$n"
         done
         ${PYTHON:-python3} "$HERE/check_kinematic_test.py" --compare "$PWD/runs/ranks_${case}_np1" \
            "$PWD/runs/ranks_${case}_np4" || FAILED+=("ranks $case: 1 and 4 ranks differ")
      done ;;
   *)
      echo "Unknown suite '$SUITE' (quick, convergence, scaling or ranks)"; exit 2 ;;
esac

if [ ${#FAILED[@]} -eq 0 ]; then
   echo "Suite '$SUITE': all runs passed"
else
   echo "Suite '$SUITE': ${#FAILED[@]} failure(s):"; printf '   %s\n' "${FAILED[@]}"; exit 1
fi
