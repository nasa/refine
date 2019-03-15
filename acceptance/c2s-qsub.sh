#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error

testname=c2s
queue=K3-route
nprocs=16
walltime=20

BUILDLOG=$PWD/../log.build
RUNLOG=$PWD/../log.runtest

cat << EOF > ${testname}.pbs
#PBS -S /bin/bash
#PBS -l select=1:ncpus=${nprocs}:mpiprocs=${nprocs}
#PBS -l walltime=00:${walltime}:00
#PBS -q ${queue}
#PBS -N ${testname}
#PBS -o ${testname}
#PBS -W umask=027
#PBS -j oe
#PBS -m n

set -o xtrace   # show commands invoked
set -o errexit  # exit on first error
set -o pipefail # ensure exit with pipe

uname -mrn

# Setup bash module environment
. /usr/local/pkgs/modules/init/bash

echo \$PBS_JOBID
cat \$PBS_NODEFILE

cd \$PBS_O_WORKDIR
source acceptance/${testname}-modules.sh

egads=\$PBS_O_WORKDIR/egads
eagds_bin=\$PBS_O_WORKDIR/egads/bin
c2s=\$PBS_O_WORKDIR/../C2S

./bootstrap
mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=\$egads \
    --with-EGADS="/u/shared/fun3d/fun3d_users/modules/ESP/114/EngSketchPad" \
    CFLAGS="-g -O2" \
    && make -j \
    && make install \
    ) \
    || exit 1

export PATH=\$PATH:\$egads_bin

cd \$c2s

dir=turbulence_modeling_resource/3D/onera-m6/geometry
(cd \$dir && ./init-grid.sh || touch \$dir/FAILED ) &

wait

find \$c2s -name FAILED

echo -e \\n\
# Build has failed if any failed cases have been reported
exit `find \$c2s -name FAILED | wc -l`

EOF

tail_file(){
  (
  while [[ ! -f $1 ]]; do
    sleep 1
  done
  tail -f -n +1 $1
  ) &
}

rm -rf ${BUILDLOG} ${RUNLOG}
(qsub -V -Wblock=true ${testname}.pbs) &
pid=$!

set +x
tail_file ${BUILDLOG}
tail_file ${RUNLOG}
set -x

# Capture test-suite error code, turn of exit-on-error
set +e 
wait ${pid}; exit_code=$?

exit ${exit_code}
