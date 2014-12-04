#!/usr/bin/env bash

bin=${HOME}/refine/strict/two

${bin}/ref_acceptance 2 ref_bamg_test.b8.ugrid

function adapt_cycle {
    proj=$1

    cp ref_bamg_test.b8.ugrid ${proj}.b8.ugrid

    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.html
    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.tec

    ${bin}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric 0.02 0.5
    ${bin}/ref_bamg ${proj}.b8.ugrid ${proj}.metric

    ${bin}/ref_metric_test ${proj}.b8.ugrid ${proj}.metric
    cp ref_metric_test_surf.tec ${proj}_metric_surf.tec
    cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
}

adapt_cycle accept-2d-bamg-00
adapt_cycle accept-2d-bamg-01
adapt_cycle accept-2d-bamg-02
adapt_cycle accept-2d-bamg-03
adapt_cycle accept-2d-bamg-04
adapt_cycle accept-2d-bamg-05
adapt_cycle accept-2d-bamg-06
adapt_cycle accept-2d-bamg-07
adapt_cycle accept-2d-bamg-08
adapt_cycle accept-2d-bamg-09

