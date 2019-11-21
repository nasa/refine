#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

serveCSM -batch hemisphere.csm

ref bootstrap hemisphere.egads

ref adapt hemisphere-vol.meshb -g hemisphere.egads -x hemicurve.meshb

