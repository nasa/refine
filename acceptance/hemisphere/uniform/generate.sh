#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

serveCSM -batch hemisphere.csm

ref_geom_test --tetgen hemisphere.egads hemicurve.meshb .5  1  15.00
