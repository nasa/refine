#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/src
    two=$1/two
else
    one=${HOME}/refine/strict/src
    two=${HOME}/refine/strict/two
fi

${two}/ref_acceptance 1 ref_adapt_test.b8.ugrid

function adapt_cycle {
    proj=$1
    field=polar-1
    
    cp ref_adapt_test.b8.ugrid ${proj}.b8.ugrid

    ${two}/ref_translate ${proj}.b8.ugrid ${proj}.html
    ${two}/ref_translate ${proj}.b8.ugrid ${proj}.tec

    ${two}/ref_acceptance -ugawg ${field} ${proj}.b8.ugrid ${proj}.metric
    
    rm ref_adapt_test.b8.ugrid
    ${two}/ref_driver -i ${proj}.b8.ugrid -m ${proj}.metric
    cp ref_driver.b8.ugrid ref_adapt_test.b8.ugrid
    
    ${two}/ref_metric_test ${proj}.b8.ugrid ${proj}.metric > ${proj}.status
    cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
}

adapt_cycle accept-3d-two-00
adapt_cycle accept-3d-two-01
adapt_cycle accept-3d-two-02
adapt_cycle accept-3d-two-03
adapt_cycle accept-3d-two-04
adapt_cycle accept-3d-two-05

cat accept-3d-two-05.status
../../../check.rb accept-3d-two-05.status 0.004 2.1



