#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

cat > faux_input <<EOF
6
1 xplane 0
2 xplane 1
3 yplane 0
4 yplane 1
5 zplane 0
6 zplane 1
EOF

if [ $# -gt 0 ] ; then
    one=$1/src
    two=$1/two
else
    one=${HOME}/refine/strict/src
    two=${HOME}/refine/strict/two
fi

${two}/ref_acceptance 1 ref_adapt_test.b8.ugrid

function adapt_cycle {
    proj=$1

    cp ref_adapt_test.b8.ugrid ${proj}.b8.ugrid

    ${two}/ref_translate ${proj}.b8.ugrid ${proj}.html
    ${two}/ref_translate ${proj}.b8.ugrid ${proj}.tec

    ${two}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric 0.001

    ${two}/ref_translate ${proj}.b8.ugrid ${proj}.fgrid
    ${one}/refine-px -g ${proj}.fgrid -m ${proj}.metric -o refine-one.fgrid
    ${two}/ref_translate refine-one.fgrid ref_adapt_test.b8.ugrid

    ${two}/ref_metric_test ${proj}.b8.ugrid ${proj}.metric > ${proj}.status
    cp ref_metric_test_surf.tec ${proj}_metric_surf.tec
    cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
    gnuplot ref_histogram_edge-ratio.gnuplot
    epstopdf ref_histogram_edge-ratio.eps
    cp ref_histogram_edge-ratio.pdf ${proj}_edge-ratio.pdf
}

adapt_cycle accept-3d-one-00
adapt_cycle accept-3d-one-01
adapt_cycle accept-3d-one-02
adapt_cycle accept-3d-one-03

cat accept-3d-one-03.status
../../../check.rb accept-3d-one-03.status 0.035 1.9

