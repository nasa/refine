#!/usr/bin/env bash

bin=${HOME}/refine/strict/two

${bin}/ref_acceptance 2 ref_bamg_test.b8.ugrid

c=0.02
k=0.5

function adapt_cycle {
    proj=$1

    cp ref_bamg_test.b8.ugrid ${proj}.b8.ugrid

    ${bin}/ref translate ${proj}.b8.ugrid ${proj}.tec

    ${bin}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric $c $k
    ${bin}/ref_bamg ${proj}.b8.ugrid ${proj}.metric

    ${bin}/ref_metric_test ${proj}.b8.ugrid ${proj}.metric
    cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
}

function gen_only {
    proj=$1

    cp ref_bamg_test.b8.ugrid ${proj}.b8.ugrid

    ${bin}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric $c $k
    ${bin}/ref_bamg ${proj}.b8.ugrid ${proj}.metric
}

gen_only accept-2d-bamg-00
gen_only accept-2d-bamg-01
gen_only accept-2d-bamg-02
gen_only accept-2d-bamg-03
gen_only accept-2d-bamg-04
gen_only accept-2d-bamg-05
gen_only accept-2d-bamg-06
gen_only accept-2d-bamg-07
gen_only accept-2d-bamg-08
adapt_cycle accept-2d-bamg-09

