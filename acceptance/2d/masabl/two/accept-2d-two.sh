#!/usr/bin/env bash

bin=${HOME}/refine/strict/src

${bin}/ref_acceptance 2 ref_adapt_test.b8.ugrid

function adapt_cycle {
    proj=$1

    cp ref_adapt_test.b8.ugrid ${proj}.b8.ugrid

    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.html
    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.tec

    ${bin}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric -masabl
    ${bin}/ref_adapt_test ${proj}.b8.ugrid ${proj}.metric | tee ${proj}.out || exit 1

    ${bin}/ref_metric_test ${proj}.b8.ugrid ${proj}.metric
    cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
}

adapt_cycle accept-2d-two-00
adapt_cycle accept-2d-two-01
adapt_cycle accept-2d-two-02
adapt_cycle accept-2d-two-03
adapt_cycle accept-2d-two-04
adapt_cycle accept-2d-two-05
adapt_cycle accept-2d-two-06
adapt_cycle accept-2d-two-07
adapt_cycle accept-2d-two-08
adapt_cycle accept-2d-two-09

