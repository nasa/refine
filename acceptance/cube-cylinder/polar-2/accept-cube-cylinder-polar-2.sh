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

    ${two}/ref adapt ${inproj}.meshb \
	  --egads ega.egads \
	  --metric ${inproj}-metric.solb \
	  -x ${outproj}.meshb -s ${sweeps} \
	  -t \
	  -f ${outproj}_stat.tec \
	  -q ${outproj}_stat.plt \
	| tee ${inproj}_refine_out
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

cat cycle02.status
../../check.rb cycle02.status 0.2 3.0

exit
