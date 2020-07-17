#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

serveCSM -batch offset.csm

ref boostrap offset.egads

ref adapt offset-vol.meshb -g offset.egads -x offset.meshb


