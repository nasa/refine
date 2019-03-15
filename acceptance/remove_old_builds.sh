#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [[ $# -ne 1 ]]; then
  echo "usage: $0 <path/to/build/prefix>"
  exit 1
fi

BUILD_PATH=$1
DAYS_OLD=14

find ${BUILD_PATH} -type d -mindepth 1 -maxdepth 1 -mtime +${DAYS_OLD} \
     -exec chmod -R u+rwX {} \;

find ${BUILD_PATH} -type d -mindepth 1 -maxdepth 1 -mtime +${DAYS_OLD} \
     -exec rm -r {} \;

