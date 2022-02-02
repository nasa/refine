#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/zoltan/one
    two=${HOME}/refine/zoltan/src
fi

mpiexec -np 1 ${two}/ref_acceptance 1 ref_adapt_test.lb8.ugrid

function adapt_cycle {
    proj=$1
    sweeps=$2
    cores=$3

    cp ref_adapt_test.lb8.ugrid ${proj}.lb8.ugrid

    mpiexec -np 1 ${two}/ref translate ${proj}.lb8.ugrid ${proj}.tec

    mpiexec -np 1 ${two}/ref_acceptance ${proj}.lb8.ugrid ${proj}.solb 0.01

    rm ref_adapt_test.lb8.ugrid
    mpiexec -np ${cores} ${two}/ref_driver -i ${proj}.lb8.ugrid -m ${proj}.solb -s ${sweeps} -x ref_adapt_test.lb8.ugrid

    mpiexec -np 1 ${two}/ref_metric_test ${proj}.lb8.ugrid ${proj}.solb > ${proj}.status
    cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
}

adapt_cycle accept-3d-linear-mpt-00 10 4
adapt_cycle accept-3d-linear-mpt-01 20 8 
adapt_cycle accept-3d-linear-mpt-02 10 8
adapt_cycle accept-3d-linear-mpt-03 10 8

cat accept-3d-linear-mpt-03.status
../../check.rb accept-3d-linear-mpt-03.status 0.3 2.0



