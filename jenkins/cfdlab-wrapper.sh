#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

git_url="git@gitlab.larc.nasa.gov:fun3d-developers/refine.git"

build_directory_root=/ssd/fun3d/jenkins

build_machine=cmb20
ssh -o StrictHostKeyChecking=no fun3d@${build_machine} true

ssh fun3d@${build_machine} <<EOF
whoami && \
cd ${build_directory_root} && \
  mkdir -p ${BUILD_TAG} && \
  cd ${BUILD_TAG} && \
    pwd && \
    git clone ${git_url} && \
    cd refine && \
      pwd && \
      git checkout '${GIT_COMMIT}'
EOF

# when merging:
#      git checkout '${gitlabSourceBranch}' && \
#      git merge '${gitlabTargetBranch}'


ssh fun3d@${build_machine} <<EOF
whoami && \
cd ${build_directory_root} && \
  cd ${BUILD_TAG} && \
    cd refine && \
      pwd && \
      git status'
EOF

ssh fun3d@${build_machine} <<EOF
whoami && \
cd ${build_directory_root} && \
  cd ${BUILD_TAG} && \
    cd refine && \
      pwd && \
      ./jenkins/cfdlab.sh
EOF

ssh -o StrictHostKeyChecking=no fun3d@${build_machine} true

scp fun3d@${build_machine}:${build_directory_root}/${BUILD_TAG}/log.\* .
scp fun3d@${build_machine}:${build_directory_root}/${BUILD_TAG}/refine-\*.tar.gz .

ssh fun3d@${build_machine} <<EOF
 whoami && \
 cd ${build_directory_root}/${BUILD_TAG}/refine && \
  ./jenkins/remove_old_builds.sh \
   ${BUILD_NUMBER} \
   "${build_directory_root}/jenkins-${JOB_NAME}"
EOF

