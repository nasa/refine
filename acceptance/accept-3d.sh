#!/usr/bin/env bash

bin=${HOME}/refine/strict/two

${bin}/ref_acceptance 1 ref_adapt_test.b8.ugrid

function adapt_cycle {
    proj=$1

    cp ref_adapt_test.b8.ugrid ${proj}.b8.ugrid

    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.html
    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.tec

    ${bin}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric 0.005
    ${bin}/ref_adapt_test ${proj}.b8.ugrid ${proj}.metric

}

adapt_cycle accept-3d-00
adapt_cycle accept-3d-01
adapt_cycle accept-3d-02
adapt_cycle accept-3d-03
adapt_cycle accept-3d-04
adapt_cycle accept-3d-05
adapt_cycle accept-3d-06

