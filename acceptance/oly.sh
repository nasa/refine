#!/usr/bin/env bash

bin=${HOME}/refine/strict/two

${bin}/ref_olympics 
${bin}/ref_translate ref_olympics.meshb  ref_adapt_test.b8.ugrid

function adapt_cycle {
    proj=$1

    cp ref_adapt_test.b8.ugrid ${proj}.b8.ugrid

    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.html
    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.meshb

    ${bin}/ref_olympics ${proj}.b8.ugrid ${proj}.metric 0.01
    ${bin}/ref_adapt_test ${proj}.b8.ugrid ${proj}.metric

}

adapt_cycle oly00
adapt_cycle oly01
adapt_cycle oly02
adapt_cycle oly03
adapt_cycle oly04
adapt_cycle oly05
