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
field="-u sin50xy"
egads="-g square.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2
    complexity=$3

    ${src}/ref_acceptance ${field} ${inproj}.meshb \
	  ${inproj}.solb

    ${src}/ref multiscale ${inproj}.meshb ${inproj}.solb \
	  ${complexity} ${inproj}-metric.solb --norm-power 1

    ${src}/ref adapt ${inproj}.meshb ${egads} -m ${inproj}-metric.solb \
	  -x ${outproj}.meshb -f ${outproj}.tec

    ${src}/ref_acceptance ${field} ${outproj}.meshb \
	  ${outproj}.solb

    ${src}/ref visualize ${outproj}.meshb ${outproj}.solb ${outproj}.plt
}

function disc_interp_error {
    proj=$1
    ${src}/ref_geom_test --enrich2 ${proj}.meshb square.egads
    mv ref_geom_enrich2.meshb ${proj}-enrich2.meshb
    ${src}/ref_acceptance ${field} ${proj}-enrich2.meshb \
	  ${proj}-enrich2.solb
    ${src}/ref_metric_test --error2 ${proj}-enrich2.meshb ${proj}-enrich2.solb
    cp ref_metric_test_dist.plt ${proj}-disc-err.plt
    ${src}/ref_metric_test --imply ${proj}.meshb ${proj}-imply.solb
    ${src}/ref_acceptance ${field} ${proj}.meshb ${proj}.solb
    ${src}/ref_metric_test --error ${proj}.meshb ${proj}.solb ${proj}-imply.solb
    cp ref_metric_test_error.plt ${proj}-cont-err.plt
}

cp square.meshb cycle00.meshb

adapt_cycle cycle00 cycle01 1000
disc_interp_error cycle00
adapt_cycle cycle01 cycle02 1000
disc_interp_error cycle01
adapt_cycle cycle02 cycle03 2000
disc_interp_error cycle02
adapt_cycle cycle03 cycle04 4000
disc_interp_error cycle04
adapt_cycle cycle04 cycle05 8000
disc_interp_error cycle05


