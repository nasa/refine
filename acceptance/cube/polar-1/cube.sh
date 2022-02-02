#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

serveCSM -batch cube.csm

ref boostrap cube.egads

ref adapt cube-vol.meshb -g cube.egads -x cube.meshb


