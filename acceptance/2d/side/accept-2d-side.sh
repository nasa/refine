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

${two}/ref_acceptance 2 ref_adapt_test.meshb

function adapt_cycle {
    proj=$1

    cp ref_adapt_test.meshb ${proj}.meshb

    ${two}/ref_translate ${proj}.meshb ${proj}.tec

    ${two}/ref_acceptance -twod side ${proj}.meshb ${proj}.solb
    ${two}/ref_driver -i ${proj}.meshb -s 1 -m ${proj}.solb -x ref_adapt_test.meshb -t | tee ${proj}.out || exit 1

    cp ref_gather_movie.tec ${proj}_movie.tec
    cp ref_gather_histo.tec ${proj}_histo.tec

    ${two}/ref_metric_test ${proj}.meshb ${proj}.solb > ${proj}.status
    cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
}

adapt_cycle accept-2d-side-00

exit

adapt_cycle accept-2d-side-01
adapt_cycle accept-2d-side-02
adapt_cycle accept-2d-side-03
adapt_cycle accept-2d-side-04

cat accept-2d-side-04.status
../../check.rb accept-2d-side-04.status 0.40 1.6

