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

field=polar-2

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3
    cores=$4

    mpiexec -np ${cores} ${two}/ref_driver -i ${inproj}.meshb -g ega.egads -m ${inproj}.metric -o ${outproj} -s ${sweeps} -t
    mv ref_gather_movie.tec ${inproj}_movie.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.metric
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.metric > ${outproj}.status
}

# ${two}/ref_geom_test ega.egads
# ${two}/ref_geom_test ega.egads ega.ugrid

${two}/ref_acceptance -ugawg ${field} ega.meshb ega.metric

adapt_cycle ega para01 2 1
adapt_cycle para01 para02 5 1
adapt_cycle para02 para03 10 4
adapt_cycle para03 para04 15 4

cat para04.status
../../../check.rb para04.status 0.1 3.5

exit
