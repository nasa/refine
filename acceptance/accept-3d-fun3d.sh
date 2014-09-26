#!/usr/bin/env bash

cat > faux_input <<EOF
6
1 xplane 0
2 xplane 1
3 yplane 0
4 yplane 1
5 zplane 0
6 zplane 1
EOF

two=${HOME}/refine/strict/two

${two}/ref_acceptance 1 ref_adapt_test.lb8.ugrid

function adapt_cycle {
    proj=$1

     cp ref_adapt_test.lb8.ugrid ${proj}.lb8.ugrid

#     ${two}/ref_translate ${proj}.fgrid ${proj}.html
    ${two}/ref_translate ${proj}.lb8.ugrid ${proj}.tec

    ${two}/ref_acceptance ${proj}.lb8.ugrid ${proj}.metric 0.001

cat > ${proj}.mapbc <<EOF
6
1 5000
2 5000
3 5000
4 5000
5 5000
6 5000
EOF

cat > fun3d.nml <<EOF
&project
project_rootname = '${proj}'
/
&raw_grid
grid_format='aflr3'
data_format='stream'
/
&adapt_metric_construction
adapt_metric_from_file='${proj}.metric'
/
&reference_physical_properties
mach_number = 1.0
reynolds_number = 1.0
/
&code_run_control
restart_read='off'
/
&adapt_mechanics
adapt_project='ref_adapt_test'
/
EOF

    nodet_mpi --adapt

}

adapt_cycle accept-3d-00
adapt_cycle accept-3d-01
adapt_cycle accept-3d-02
adapt_cycle accept-3d-03
adapt_cycle accept-3d-04
adapt_cycle accept-3d-05
adapt_cycle accept-3d-06

