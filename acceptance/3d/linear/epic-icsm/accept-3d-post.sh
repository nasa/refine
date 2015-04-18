#!/usr/bin/env bash

bin=${HOME}/refine/strict/two

proj=epic_L2_qShape_InsertCollapseSwapMove

${bin}/ref_translate ${proj}.b8.ugrid accept-3d.tec

${bin}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric 0.001

${bin}/ref_metric_test ${proj}.b8.ugrid ${proj}.metric > ${proj}.status
cp ref_metric_test_surf.tec ${proj}_metric_surf.tec
cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
gnuplot ref_histogram_edge-ratio.gnuplot
epstopdf ref_histogram_edge-ratio.eps
cp ref_histogram_edge-ratio.pdf ${proj}_edge-ratio.pdf
