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
field1="-u combo1"
field2="-u combo2"
egads="-g square.egads"

function adapt_cycle {
    inproj=$1
    outproj=$2
    complexity=$3

    ${src}/ref_acceptance ${field1} ${inproj}.meshb \
	  ${inproj}-1.solb
    ${src}/ref_acceptance ${field2} ${inproj}.meshb \
 	  ${inproj}-2.solb

    ${src}/ref metric \
	  ${inproj}.meshb ${inproj}-1.solb \
	  ${complexity} ${inproj}-metric1.solb
    ${src}/ref metric \
	  ${inproj}.meshb ${inproj}-2.solb \
	  ${complexity} ${inproj}-metric2.solb

    ${src}/ref_metric_test --error \
	  ${inproj}.meshb ${inproj}-1.solb ${inproj}-metric1.solb
    cp ref_metric_test_error.plt ${inproj}-metric1-error1.plt
    ${src}/ref_metric_test --error \
	  ${inproj}.meshb ${inproj}-1.solb ${inproj}-metric2.solb
    cp ref_metric_test_error.plt ${inproj}-metric2-error1.plt
    ${src}/ref_metric_test --error \
	  ${inproj}.meshb ${inproj}-2.solb ${inproj}-metric1.solb
    cp ref_metric_test_error.plt ${inproj}-metric1-error2.plt
    ${src}/ref_metric_test --error \
	  ${inproj}.meshb ${inproj}-2.solb ${inproj}-metric2.solb
    cp ref_metric_test_error.plt ${inproj}-metric2-error2.plt

    ${src}/ref_metric_test --combine \
	  ${inproj}.meshb ${inproj}-1.solb ${inproj}-2.solb \
	  2 -1 ${complexity} 0.5 ${inproj}-metric.solb

    ${src}/ref_metric_test --error \
	  ${inproj}.meshb ${inproj}-1.solb ${inproj}-metric.solb
    cp ref_metric_test_error.plt ${inproj}-metric-error1.plt
    ${src}/ref_metric_test --error \
	  ${inproj}.meshb ${inproj}-2.solb ${inproj}-metric.solb
    cp ref_metric_test_error.plt ${inproj}-metric-error2.plt

    ${src}/ref adapt ${inproj}.meshb ${egads} -m ${inproj}-metric.solb \
	  -x ${outproj}.meshb -f ${outproj}.tec

    ${src}/ref_acceptance ${field1} ${outproj}.meshb \
	  ${outproj}-1.solb
    ${src}/ref_acceptance ${field2} ${outproj}.meshb \
 	  ${outproj}-2.solb

    ${src}/ref visualize ${outproj}.meshb ${outproj}-1.solb ${outproj}-1.plt
    ${src}/ref visualize ${outproj}.meshb ${outproj}-2.solb ${outproj}-2.plt
}

cp square.meshb cycle00.meshb

adapt_cycle cycle00 cycle01 1000
adapt_cycle cycle01 cycle02 1000
adapt_cycle cycle02 cycle03 2000
adapt_cycle cycle03 cycle04 4000
#adapt_cycle cycle04 cycle05 8000


