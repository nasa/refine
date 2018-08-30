#!/usr/bin/env bash

bin=${HOME}/refine/strict/src

${bin}/ref_acceptance 2 ref_adapt_test.b8.ugrid

function adapt_cycle {
    proj=$1

    cp ref_adapt_test.b8.ugrid ${proj}.b8.ugrid

    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.html
    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.tec

    ${bin}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric -masabl
    ${bin}/ref_driver -i ${proj}.b8.ugrid -m ${proj}.metric -o ref_adapt_test -t | tee ${proj}.out || exit 1

    cp ref_gather_movie.tec ${proj}_movie.tec
    cp ref_gather_histo.tec ${proj}_histo.tec

    ${bin}/ref_metric_test ${proj}.b8.ugrid ${proj}.metric > ${proj}.status
    cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
}

adapt_cycle accept-2d-two-00
adapt_cycle accept-2d-two-01
adapt_cycle accept-2d-two-02
adapt_cycle accept-2d-two-03
adapt_cycle accept-2d-two-04
adapt_cycle accept-2d-two-05

cat accept-2d-two-05.status
../../../check.rb accept-2d-two-05.status 0.45 1.6


