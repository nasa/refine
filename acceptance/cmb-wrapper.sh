#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

build_directory_root=/scratch/fun3d/gitlab-ci

build_machine=cmb20
ssh -o StrictHostKeyChecking=no fun3d@${build_machine} true

BUILD_TAG="${CI_JOB_NAME}-${CI_JOB_ID}"
set +u
if [[ -z "${CI_MERGE_REQUEST_SOURCE_BRANCH_NAME}" || -z "{CI_MERGE_REQUEST_TARGET_BRANCH_NAME}" ]]; then
    checkout_cmd="git checkout ${CI_COMMIT_SHA}"
else
    checkout_cmd="git checkout ${CI_MERGE_REQUEST_SOURCE_BRANCH_NAME} && git merge ${CI_MERGE_REQUEST_TARGET_BRANCH_NAME}"
fi
set -u

ssh -o LogLevel=error fun3d@${build_machine} <<EOF
whoami && \
mkdir -p ${build_directory_root} && \
cd ${build_directory_root} && \
  mkdir -p ${BUILD_TAG} && \
  cd ${BUILD_TAG} && \
    pwd && \
    git clone ${CI_REPOSITORY_URL} && \
    cd refine && \
      pwd && \
      ${checkout_cmd} && \
    ./acceptance/cmb.sh
EOF

scp -o LogLevel=error fun3d@${build_machine}:${build_directory_root}/${BUILD_TAG}/log.\* .
scp -o LogLevel=error fun3d@${build_machine}:${build_directory_root}/${BUILD_TAG}/refine-\*.tar.gz .

trap "cat cleanup.log" EXIT
ssh -o LogLevel=error fun3d@${build_machine} > cleanup.log 2>&1 <<EOF
whoami && \
uname -n && \
cd ${build_directory_root}/${BUILD_TAG}/refine && \
 ./acceptance/remove_old_builds.sh \
  "${build_directory_root}"
EOF
trap - EXIT

