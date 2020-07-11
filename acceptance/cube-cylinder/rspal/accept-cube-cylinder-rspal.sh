#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

field=polar-2

function adapt_cycle {
    inproj=$1
    outproj=$2
    sweeps=$3

    egads="-g cube-cylinder.egads"

    ${src}/ref_acceptance -ugawg ${field} ${inproj}.meshb ${inproj}-metric.solb
    ${src}/ref adapt ${inproj}.meshb ${egads} -m ${inproj}-metric.solb -x ${outproj}.meshb -s ${sweeps} -t -f ${outproj}_stat.tec | tee ${outproj}_refine_out
    cp ref_gather_movie.tec ${outproj}_movie.tec
    cp ref_gather_histo.tec ${outproj}_histo.tec
    ${src}/ref_acceptance -ugawg ${field} ${outproj}.meshb ${outproj}-metric.solb
    ${src}/ref_metric_test ${outproj}.meshb ${outproj}-metric.solb > ${outproj}.status
}

adapt_cycle cube-cylinder00 cube-cylinder01 2
adapt_cycle cube-cylinder01 cube-cylinder02 5
adapt_cycle cube-cylinder02 cube-cylinder03 5
adapt_cycle cube-cylinder03 cube-cylinder04 5

cat cube-cylinder04.status
../../check.rb cube-cylinder04.status 0.2 3.0

exit
