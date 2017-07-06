#!/usr/bin/env bash

bin=${HOME}/refine/strict/two

proj=epic_L2_qShape_InsertCollapseSwap

${bin}/ref_translate ${proj}.b8.ugrid accept-3d.tec

${bin}/ref_acceptance ${proj}.b8.ugrid ${proj}.metric 0.001

${bin}/ref_metric_test ${proj}.b8.ugrid ${proj}.metric > ${proj}.status
cp ref_metric_test_s00_n1_p0_ellipse.tec ${proj}_metric_ellipse.tec
