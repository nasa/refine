#!/usr/bin/env bash

bin=${HOME}/refine/strict/two

${bin}/ref_acceptance 2 ref_bamg_test.b8.ugrid

function adapt_cycle {
    proj=$1

    cp ref_bamg_test.b8.ugrid ${proj}.b8.ugrid

    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.html
    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.tec

    ${bin}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric 0.0001
    ${bin}/ref_bamg ${proj}.b8.ugrid ${proj}.metric

}

adapt_cycle accept-bamg-00
adapt_cycle accept-bamg-01
adapt_cycle accept-bamg-02
adapt_cycle accept-bamg-03
adapt_cycle accept-bamg-04
adapt_cycle accept-bamg-05
adapt_cycle accept-bamg-06
adapt_cycle accept-bamg-07
adapt_cycle accept-bamg-08
adapt_cycle accept-bamg-09

