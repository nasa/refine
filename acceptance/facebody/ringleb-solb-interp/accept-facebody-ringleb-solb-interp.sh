#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

tecplot=-t
field="-q ringleb"
interp="-u u5"
egads="--egads ringleb.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2
    complexity=$3

    ${src}/ref_acceptance ${field} ${inproj}.meshb \
	  ${inproj}_volume.solb
    ${src}/ref_gather_test ${inproj}.meshb \
	  ${inproj}_volume.solb ${inproj}_volume.tec

    ${src}/ref_acceptance ${interp} ${inproj}.meshb \
	  ${inproj}-interp.solb
    ${src}/ref loop ${inproj} ${outproj} ${complexity} \
	  --interpolant ${inproj}-interp.solb \
	  ${egads} -s 5 > ${inproj}-loop.txt

    ${src}/ref_acceptance ${field} ${outproj}.meshb \
	  ${outproj}_volume.solb
    ${src}/ref_gather_test ${outproj}.meshb \
	  ${outproj}_volume.solb ${outproj}_volume.tec

}

serveCSM -batch ringleb.csm > ringleb-servecsm.txt
${src}/ref bootstrap ringleb.egads > ringleb-bootstrap.txt
mv ringleb-vol.meshb cycle00.meshb

adapt_cycle cycle00 cycle01 1000
adapt_cycle cycle01 cycle02 1000
adapt_cycle cycle02 cycle03 2000
adapt_cycle cycle03 cycle04 2000
adapt_cycle cycle04 cycle05 4000
adapt_cycle cycle05 cycle06 4000

../../check.rb cycle05-loop.txt 0.2 2.0

