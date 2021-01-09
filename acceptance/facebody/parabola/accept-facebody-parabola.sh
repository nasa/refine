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
field="parabola"
egads="-g square.egads"

function func_cycle {
    inproj=$1
    outproj=$2

    ${src}/ref_acceptance -u ${field} ${inproj}.meshb \
	  ${inproj}.solb

    ${src}/ref multiscale ${inproj}.meshb ${inproj}.solb \
	  1000 ${inproj}-metric.solb

    ${src}/ref adapt ${inproj}.meshb ${egads} -m ${inproj}-metric.solb \
	  -x ${outproj}.meshb -f ${outproj}.tec

    ${src}/ref_acceptance -u ${field} ${outproj}.meshb \
	  ${outproj}.solb
}

function grad_cycle {
    inproj=$1
    outproj=$2

    ${src}/ref_acceptance -g ${field} ${inproj}.meshb \
	  ${inproj}.solb

    ${src}/ref_metric_test --multigrad ${inproj}.meshb ${inproj}.solb \
	  2 -1 1000 ${inproj}-metric.solb

    ${src}/ref adapt ${inproj}.meshb ${egads} -m ${inproj}-metric.solb \
	  -x ${outproj}.meshb -f ${outproj}.tec

    ${src}/ref_acceptance -g ${field} ${outproj}.meshb \
	  ${outproj}.solb
}

serveCSM -batch square.csm
${src}/ref bootstrap square.egads
mv square-vol.meshb func00.meshb

func_cycle func00 func01
func_cycle func01 func02
func_cycle func02 func03
func_cycle func03 func04

cp func00.meshb grad00.meshb

grad_cycle grad00 grad01
grad_cycle grad01 grad02
grad_cycle grad02 grad03
grad_cycle grad03 grad04


