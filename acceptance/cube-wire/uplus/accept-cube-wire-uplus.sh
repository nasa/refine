#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

egads="-g cube-wire.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2
    complexity=$3
    sweeps=5

    ${src}/ref \
	  adapt ${inproj}.meshb \
	  ${egads} \
	  --spalding 0.001 ${complexity} \
	  --viscous-tags 7 \
	  -x ${outproj}.meshb -s ${sweeps} \
	  -t -f ${outproj}_stat.tec | tee ${outproj}_refine_out
    cp ref_gather_movie.tec ${outproj}_movie.tec
    cp ref_gather_histo.tec ${outproj}_histo.tec
}

rm -f cube-wire.egads
serveCSM -batch cube-wire.csm
rm -f cube-wire-vol.meshb
${src}/ref bootstrap cube-wire.egads
mv cube-wire-vol.meshb cube-wire00.meshb

adapt_cycle cube-wire00 cube-wire01 1000
adapt_cycle cube-wire01 cube-wire02 1000

exit
