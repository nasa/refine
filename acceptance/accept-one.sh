#!/usr/bin/env bash

two=${HOME}/refine/strict/two
one=${HOME}/refine/strict/src

${two}/ref_acceptance 1 ref_adapt_test.b8.ugrid
${two}/ref_translate ref_adapt_test.b8.ugrid refine-one.fgrid

function adapt_cycle {
    proj=$1

    cp refine-one.fgrid ${proj}.fgrid

#     ${two}/ref_translate ${proj}.fgrid ${proj}.html
    ${two}/ref_translate ${proj}.fgrid ${proj}.tec

    ${two}/ref_acceptance ${proj}.fgrid ${proj}.metric 0.001
    ${one}/refine-px -g ${proj}.fgrid -m ${proj}.metric -o refine-one.fgrid

}

adapt_cycle accept-3d-00
adapt_cycle accept-3d-01
adapt_cycle accept-3d-02
adapt_cycle accept-3d-03
adapt_cycle accept-3d-04
adapt_cycle accept-3d-05
adapt_cycle accept-3d-06

