#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

serveCSM -batch sphcyl.csm

ref boostrap sphcyl.egads --blend sphcyl-blend.meshb

ref adapt sphcyl-vol.meshb -g sphcyl.egads --blend sphcyl-blend.meshb -x sphcyl.meshb

