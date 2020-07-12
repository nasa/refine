#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

field="-u cyluplus"
egads="-g cube-cylinder.egads"

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

adapt_cycle cube-cylinder00 cube-cylinder01 500
adapt_cycle cube-cylinder01 cube-cylinder02 1000
adapt_cycle cube-cylinder02 cube-cylinder03 2000
adapt_cycle cube-cylinder03 cube-cylinder04 4000
adapt_cycle cube-cylinder04 cube-cylinder05 8000
adapt_cycle cube-cylinder05 cube-cylinder06 16000
adapt_cycle cube-cylinder06 cube-cylinder07 32000

exit
