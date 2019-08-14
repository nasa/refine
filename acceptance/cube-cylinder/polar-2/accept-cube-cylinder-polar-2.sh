#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/egads/one
    two=${HOME}/refine/egads/src
fi

field=polar-2

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    ${two}/ref_driver -i ${inproj}.meshb -g ega.egads -m ${inproj}-metric.solb -o ${outproj} -s ${sweeps} -t | tee ${inproj}_refine_out
    cp ref_gather_movie.tec ${inproj}_movie.tec
    cp ref_gather_histo.tec ${inproj}_histo.tec
    ${two}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}-metric.solb
    ${two}/ref_metric_test ${outproj}.meshb ${outproj}-metric.solb > ${outproj}.status
}

# ${two}/ref_geom_test ega.egads
# ${two}/ref_geom_test ega.egads ega.ugrid

${two}/ref_acceptance -ugawg ${field} ega.meshb ega-metric.solb

adapt_cycle ega cycle01 2
adapt_cycle cycle01 cycle02 15
adapt_cycle cycle02 cycle03 15
adapt_cycle cycle03 cycle04 15

cat cycle04.status
../../check.rb cycle04.status 0.2 3.0

exit
