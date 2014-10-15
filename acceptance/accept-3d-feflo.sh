#!/usr/bin/env bash

bin=${HOME}/refine/strict/two

${bin}/ref_acceptance 1 ref_adapt_test.b8.ugrid

function adapt_cycle {
    proj=$1

    cp ref_adapt_test.b8.ugrid ${proj}.b8.ugrid

    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.html
    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.tec

    ${bin}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric 0.001

    ${bin}/ref_translate ${proj}.b8.ugrid ${proj}.meshb
    ~/refine/util/metric2sol.sh ${proj}.metric ${proj}.sol

    ~/inria-file-tools/mac-bin/feflo.a -in ${proj}.meshb -met ${proj}.sol -out ref_adapt_test.meshb

    ${bin}/ref_translate ref_adapt_test.meshb ref_adapt_test.b8.ugrid
    ${bin}/ref_metric_test ${proj}.b8.ugrid ${proj}.metric
    cp ref_metric_test_surf.tec ${proj}_metric_surf.tec
    cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
    gnuplot ref_histogram_edge-ratio.gnuplot
    epstopdf ref_histogram_edge-ratio.eps
    cp ref_histogram_edge-ratio.pdf ${proj}_edge-ratio.pdf
}

adapt_cycle accept-3d-feflo-00
adapt_cycle accept-3d-feflo-01
adapt_cycle accept-3d-feflo-02
adapt_cycle accept-3d-feflo-03
adapt_cycle accept-3d-feflo-04
adapt_cycle accept-3d-feflo-05
adapt_cycle accept-3d-feflo-06

