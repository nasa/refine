#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

field="-u sphuplus"
egads="-g cube-sphere.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2
    complexity=$3
    sweeps=10

    ${src}/ref_acceptance ${field} ${inproj}.meshb ${inproj}-uplus.solb

    ${src}/ref_gather_test ${inproj}.meshb \
	  ${inproj}-uplus.solb ${inproj}-uplus.tec
    
    ${src}/ref multiscale ${inproj}.meshb ${inproj}-uplus.solb \
	  ${complexity} ${inproj}-metric.solb \
	  --gradation 10 | tee ${inproj}-multi.txt
    
    ${src}/ref \
	  adapt ${inproj}.meshb \
	  ${egads} \
	  -m ${inproj}-metric.solb \
	  -x ${outproj}.meshb -s ${sweeps} \
	  -t -f ${outproj}_stat.tec | tee ${outproj}_refine_out
    cp ref_gather_movie.tec ${outproj}_movie.tec
    cp ref_gather_histo.tec ${outproj}_histo.tec

    ${src}/ref_acceptance ${field} ${outproj}.meshb ${outproj}-uplus.solb

    ${src}/ref_gather_test ${outproj}.meshb \
	  ${outproj}-uplus.solb ${outproj}-uplus.tec
}

adapt_cycle cube-sphere00 cube-sphere01 500
adapt_cycle cube-sphere01 cube-sphere02 1000
adapt_cycle cube-sphere02 cube-sphere03 2000
adapt_cycle cube-sphere03 cube-sphere04 4000
adapt_cycle cube-sphere04 cube-sphere05 8000
adapt_cycle cube-sphere05 cube-sphere06 16000
adapt_cycle cube-sphere06 cube-sphere07 32000

exit
