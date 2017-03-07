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

${two}/ref_acceptance 2 ref_adapt_test.b8.ugrid

function adapt_cycle {
    proj=$1

    cp ref_adapt_test.b8.ugrid ${proj}.b8.ugrid

    ${two}/ref_translate ${proj}.b8.ugrid ${proj}.html
    ${two}/ref_translate ${proj}.b8.ugrid ${proj}.tec

    ${two}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric 0.0001
    ${two}/ref_adapt_test ${proj}.b8.ugrid ${proj}.metric | tee ${proj}.out || exit 1

    ${two}/ref_metric_test ${proj}.b8.ugrid ${proj}.metric > ${proj}.status
    cp ref_metric_test_surf.tec ${proj}_metric_surf.tec
    cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
    gnuplot ref_histogram_edge-ratio.gnuplot
    epstopdf ref_histogram_edge-ratio.eps
    cp ref_histogram_edge-ratio.pdf ${proj}_edge-ratio.pdf
}

adapt_cycle accept-2d-two-00
adapt_cycle accept-2d-two-01
adapt_cycle accept-2d-two-02
adapt_cycle accept-2d-two-03
adapt_cycle accept-2d-two-04
adapt_cycle accept-2d-two-05
adapt_cycle accept-2d-two-06

