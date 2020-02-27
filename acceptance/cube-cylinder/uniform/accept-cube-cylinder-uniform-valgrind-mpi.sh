#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/zoltan/one
    two=${HOME}/refine/zoltan/src
fi

nproc=4
geomfile=ega.egads

# ${two}/ref_geom_test ${geomfile}
# ${two}/ref_geom_test ${geomfile} ega.meshb

${two}/ref_acceptance ega.meshb ega.metric 0.1

mpiexec -np ${nproc} \
valgrind --quiet  --error-exitcode=1 --leak-check=full \
--gen-suppressions=all \
--suppressions=../../../misc/valgrind_suppressions_occ \
--suppressions=../../../misc/valgrind_suppressions_intel_17 \
--suppressions=../../../misc/valgrind_suppressions_openmpi \
--suppressions=../../../misc/valgrind_suppressions_zoltan \
         ${two}/ref_driver -i ega.meshb -g ${geomfile} -m ega.metric -o ref_driver1 -s 2
