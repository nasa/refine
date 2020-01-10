#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error
set -x # echo

uname -n

pwd # echo path

testname=boot

# Use this script to find quickest-to-start queue when K is backed up
queue=$(./acceptance/k-best-queue.py)

if [[ "$queue" == "K4-"* ]]; then
  nprocs=40
  walltime=50
elif [[ "$queue" == "K3-"* ]]; then
  nprocs=16
  walltime=50
elif [[ "$queue" == "K2-"* ]]; then
  nprocs=12
  walltime=50
elif [[ "$queue" == "K2a-"* ]]; then
  nprocs=12
  walltime=50
else
  echo "unknown queue requested: $queue"
  exit 1
fi


BUILDLOG=$PWD/../log-pbs-out.txt

cat << EOF > ${testname}.pbs
#PBS -S /bin/bash -q ${queue} -l select=1:ncpus=${nprocs}:mpiprocs=${nprocs} -l walltime=00:${walltime}:00
#PBS -N ${testname}
#PBS -o ${testname}
#PBS -W umask=027
#PBS -m n

set -o xtrace   # show commands invoked
set -o errexit  # exit on first error
set -o pipefail # ensure exit with pipe

uname -mrn

cd \$PBS_O_WORKDIR

./acceptance/boot.sh

EOF

trap "sleep 5; cat ${testname}.* ${testname}; pwd" EXIT
qsub -V -Wblock=true ${testname}.pbs
trap - EXIT
