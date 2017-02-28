#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

# ./remove_old_builds.sh ${BUILD_NUMBER} "path/to/build"

BUILD_NUMBER=$1
BUILD_PATH=$2

TARGET=${BUILD_NUMBER}
let TARGET=TARGET-20 || true # ((0)) has a non zero exit code
set +x
while [ ${TARGET} -gt 0 ]; do
    BUILD=${BUILD_PATH}-${TARGET}
    ((TARGET=TARGET-1)) || true # ((0)) has a non zero exit code
    if [ ! -d ${BUILD} ]; then 
      continue
    fi
    chmod -R u+rwx ${BUILD} 2> /dev/null
    rm -rf ${BUILD} 2> /dev/null
done
set -x

