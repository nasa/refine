
#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    src=$1/src
else
    src=${HOME}/refine/egads/src
fi

tecplot=""
field="-u dist"
egads="-g square-circle.egads"

cp square-circle.meshb cycle00.meshb

function adapt_cycle {
    inproj=$1
    outproj=$2

    ${src}/ref_acceptance ${field} ${inproj}.meshb \
	  ${inproj}.solb

    ${src}/ref_gather_test ${inproj}.meshb \
	  ${inproj}.solb ${inproj}-dist.tec

    ${src}/ref multiscale ${inproj}.meshb ${inproj}.solb \
	 100 ${inproj}-metric.solb

    ${src}/ref adapt ${inproj}.meshb ${egads} -m ${inproj}-metric.solb \
	  -x ${outproj}.meshb -f ${outproj}.tec

    ${src}/ref_acceptance ${field} ${outproj}.meshb \
	  ${outproj}.solb
}

adapt_cycle cycle00 cycle01
adapt_cycle cycle01 cycle02
adapt_cycle cycle02 cycle03
adapt_cycle cycle03 cycle04
adapt_cycle cycle04 cycle05
adapt_cycle cycle05 cycle06
adapt_cycle cycle06 cycle07
adapt_cycle cycle07 cycle08
adapt_cycle cycle08 cycle09
adapt_cycle cycle09 cycle10


