#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/strict/one
    two=${HOME}/refine/strict/src
fi

${two}/ref_acceptance 2 ref_adapt_test.b8.ugrid

function adapt_cycle {
    proj=$1

    cp ref_adapt_test.b8.ugrid ${proj}.b8.ugrid

    ${two}/ref_translate ${proj}.b8.ugrid ${proj}.tec

    ${two}/ref_acceptance -twod radial-1 ${proj}.b8.ugrid ${proj}.solb
    ${two}/ref_driver -i ${proj}.b8.ugrid -m ${proj}.solb -x ref_adapt_test.b8.ugrid -t | tee ${proj}.out || exit 1

    cp ref_gather_movie.tec ${proj}_movie.tec
    cp ref_gather_histo.tec ${proj}_histo.tec

    ${two}/ref_metric_test ${proj}.b8.ugrid ${proj}.solb > ${proj}.status
    cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
}

adapt_cycle accept-2d-radial-1-00
adapt_cycle accept-2d-radial-1-01
adapt_cycle accept-2d-radial-1-02
adapt_cycle accept-2d-radial-1-03
adapt_cycle accept-2d-radial-1-04

cat accept-2d-radial-1-04.status
../../check.rb accept-2d-radial-1-04.status 0.40 1.6

