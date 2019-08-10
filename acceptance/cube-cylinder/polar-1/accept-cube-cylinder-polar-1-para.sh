#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    two=$1/src
else
    two=${HOME}/refine/zoltan/src
fi

field=polar-1

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    mpiexec -np 4 ${two}/ref_driver -i ${inproj}.meshb -g ega.egads -m ${inproj}.solb -o ${outproj} -s ${sweeps} -t
    cp ref_gather_movie.tec ${inproj}_movie.tec
    cp ref_gather_histo.tec ${inproj}_histo.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}.solb
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}.solb > ${outproj}.status

}

# ${two}/ref_geom_test ega.egads ega.meshb

${two}/ref_acceptance -ugawg ${field} ega.meshb ega.solb

adapt_cycle ega cycle01 2
adapt_cycle cycle01 cycle02 2
adapt_cycle cycle02 cycle03 15
adapt_cycle cycle03 cycle04 15
adapt_cycle cycle04 cycle05 15

cat cycle04.status
../../check.rb cycle04.status 0.060 3.0

