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

    mpiexec -np ${cores} ${two}/ref adapt ${inproj}.meshb -g ega.egads -m ${inproj}-metric.solb -x ${outproj}.meshb -s ${sweeps} -t
    cp ref_gather_movie.tec ${inproj}_movie.tec
    cp ref_gather_histo.tec ${inproj}_histo.tec
    mpiexec -np 1 ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}-metric.solb
    mpiexec -np 1 ${two}/ref_metric_test ${outproj}.meshb ${outproj}-metric.solb > ${outproj}.status
}

# ${two}/ref_geom_test ega.egads
# ${two}/ref_geom_test ega.egads ega.ugrid

mpiexec -np 1 ${two}/ref_acceptance -ugawg ${field} ega.meshb ega-metric.solb

adapt_cycle ega mpt01 2 1
adapt_cycle mpt01 mpt02 15 4
adapt_cycle mpt02 mpt03 15 8
adapt_cycle mpt03 mpt04 15 8

cat mpt04.status
../../check.rb mpt04.status 0.20 3.0

exit
