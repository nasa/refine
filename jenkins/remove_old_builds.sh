#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [[ $# -ne 2 ]]; then
  echo "usage: $0 <build number> <path/to/build/prefix>"
  exit 1
fi

BUILD_NUMBER=$1
BUILD_PATH=$2

BUILDS_TO_KEEP=10
TARGET=${BUILD_NUMBER}
let TARGET=TARGET-${BUILDS_TO_KEEP}

for BUILD in ${BUILD_PATH}*; do
  if [[ ${BUILD#${BUILD_PATH}-} -lt $TARGET ]]; then
    chmod -R u+rwx ${BUILD}
    rm -rf ${BUILD}
  fi
done
