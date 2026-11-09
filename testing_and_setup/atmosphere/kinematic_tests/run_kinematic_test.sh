#!/bin/bash
#
# Run one MPAS-A kinematic transport test: init_atmosphere (case 14), atmosphere_model in
# kinematic mode, then check_kinematic_test.py. See README.md.
#
# Usage: run_kinematic_test.sh -c <case> -m <mesh> -d <mesh_dir> [-b <exe_dir>] [-n <ranks>]
#                              [-t <dt>] [-k <rk_order>] [-L <on|off>] [-l <nlev>] [-N <ntracers>]
#                              [-r <duration>]
#                              [-i <output_interval>] [-o <run_dir>]
#
set -e

HERE=$(cd "$(dirname "$0")" && pwd)
MPAS_ROOT=$(cd "$HERE/../../.." && pwd)
MPIRUN=${MPAS_MPIRUN:-mpirun -n}
PYTHON=${PYTHON:-python3}

CASE=; MESH=; MESH_DIR=; EXE_DIR=$MPAS_ROOT; NRANKS=4; DT=; RKORDER=2; LIMITER=on; NLEV=5; NTRACERS=5
DURATION=12_00:00:00; OUTPUT_INTERVAL=3_00:00:00; RUN_DIR=

usage() { sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }

while getopts "c:m:d:b:n:t:k:L:l:N:r:i:o:h" opt; do
   case $opt in
      c) CASE=$OPTARG ;;           m) MESH=$OPTARG ;;          d) MESH_DIR=$OPTARG ;;
      b) EXE_DIR=$OPTARG ;;        n) NRANKS=$OPTARG ;;        t) DT=$OPTARG ;;
      k) RKORDER=$OPTARG ;;        L) LIMITER=$OPTARG ;;
      l) NLEV=$OPTARG ;;           N) NTRACERS=$OPTARG ;;      r) DURATION=$OPTARG ;;
      i) OUTPUT_INTERVAL=$OPTARG ;; o) RUN_DIR=$OPTARG ;;      *) usage ;;
   esac
done
[ -n "$CASE" ] && [ -n "$MESH" ] && [ -n "$MESH_DIR" ] || usage
# config_time_integration_order: 2 (MPAS default; stages dt/2, dt/2, dt) or 3 (dt/3, dt/2, dt)
[ "$RKORDER" = 2 ] || [ "$RKORDER" = 3 ] || { echo "-k must be 2 or 3"; exit 2; }
# config_monotonic: FCT limiter on the last RK stage
case $LIMITER in on) MONOTONIC=true ;; off) MONOTONIC=false ;; *) echo "-L must be on or off"; exit 2 ;; esac

TEMPLATES=$HERE/$CASE
[ -d "$TEMPLATES" ] || { echo "Unknown case '$CASE' (no directory $TEMPLATES)"; exit 2; }
MESH_DIR=$(cd "$MESH_DIR" && pwd)
EXE_DIR=$(cd "$EXE_DIR" && pwd)
[ -f "$MESH_DIR/$MESH.grid.nc" ] || { echo "Missing $MESH_DIR/$MESH.grid.nc"; exit 2; }
for exe in init_atmosphere_model atmosphere_model; do
   [ -x "$EXE_DIR/$exe" ] || { echo "Missing executable $EXE_DIR/$exe"; exit 2; }
done

# Default time step: Courant number about 0.35 for the 90 m/s peak wind
if [ -z "$DT" ]; then
   case $MESH in
      x1.10242)  DT=900.0 ;;
      x1.40962)  DT=450.0 ;;
      x1.163842) DT=225.0 ;;
      *) echo "No default time step for mesh $MESH; give one with -t"; exit 2 ;;
   esac
fi

RUN_DIR=${RUN_DIR:-$PWD/runs/${CASE}_${MESH}_N${NTRACERS}_rk${RKORDER}_lim${LIMITER}}
mkdir -p "$RUN_DIR"
RUN_DIR=$(cd "$RUN_DIR" && pwd)
cd "$RUN_DIR"
find . -maxdepth 1 \( -name 'history.*.nc' -o -name 'log.*' -o -name '*.init.nc' -o -name 'results.json' \) -delete

# Mesh and partition
ln -sf "$MESH_DIR/$MESH.grid.nc" .
if [ "$NRANKS" -gt 1 ]; then
   PART=$MESH.graph.info.part.$NRANKS
   if [ -f "$MESH_DIR/$PART" ]; then
      ln -sf "$MESH_DIR/$PART" .
   elif command -v gpmetis > /dev/null; then
      cp "$MESH_DIR/$MESH.graph.info" . && gpmetis -minconn -contig -niter=200 $MESH.graph.info $NRANKS > gpmetis.log
   else
      echo "No partition file $PART and gpmetis is not available"; exit 2
   fi
fi

# Namelists and streams from the case templates
for f in namelist.init_atmosphere streams.init_atmosphere namelist.atmosphere streams.atmosphere; do
   sed -e "s/@MESH@/$MESH/g" -e "s/@NLEV@/$NLEV/g" -e "s/@NTRACERS@/$NTRACERS/g" \
       -e "s/@DT@/$DT/g" -e "s/@RKORDER@/$RKORDER/g" -e "s/@MONOTONIC@/$MONOTONIC/g" -e "s/@DURATION@/$DURATION/g" -e "s/@OUTPUT_INTERVAL@/$OUTPUT_INTERVAL/g" \
       -e "s/@CASE@/$CASE/g" "$TEMPLATES/$f" > $f
done
cp "$TEMPLATES/stream_list.atmosphere.output" .

echo "== $CASE on $MESH: dt = $DT s, RK order $RKORDER, limiter $LIMITER, $NLEV levels, $NTRACERS tracers, $DURATION, $NRANKS ranks"
echo "   run directory: $RUN_DIR"

$MPIRUN $NRANKS "$EXE_DIR/init_atmosphere_model" > init.stdout 2>&1 || { echo "init_atmosphere_model failed; see $RUN_DIR/init.stdout and log.init_atmosphere.*"; exit 1; }
mkdir -p init_logs && mv log.init_atmosphere.* init_logs/

START=$(date +%s)
$MPIRUN $NRANKS "$EXE_DIR/atmosphere_model" > atm.stdout 2>&1 || { echo "atmosphere_model failed; see $RUN_DIR/atm.stdout and log.atmosphere.*"; exit 1; }
echo "   atmosphere_model: $(( $(date +%s) - START )) s"

$PYTHON "$HERE/check_kinematic_test.py" "$RUN_DIR"
