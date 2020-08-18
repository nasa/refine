#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

serveCSM -batch cylcyl.csm

ref boostrap cylcyl.egads --blend cylcyl-blend.meshb

ref adapt cylcyl-vol.meshb -g cylcyl.egads --blend cylcyl-blend.meshb -x cylcyl.meshb

