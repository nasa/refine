#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

git_url="git@gitlab.larc.nasa.gov:fun3d-developers/refine.git"

build_directory_root=/work13/fun3d/jenkins

build_machine=cypher-work13
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
      git checkout '${gitlabSourceBranch}' && \
      git merge '${gitlabTargetBranch}'
EOF

build_machine="i17n1"
ssh -o StrictHostKeyChecking=no fun3d@${build_machine} true

ssh fun3d@${build_machine} <<EOF
whoami && \
cd ${build_directory_root} && \
  cd ${BUILD_TAG} && \
    cd refine && \
      pwd && \
      ./jenkins/cfdlab.sh
EOF

build_machine=cypher-work13
ssh -o StrictHostKeyChecking=no fun3d@${build_machine} true

ssh fun3d@${build_machine} <<EOF
 whoami && \
 cd ${build_directory_root}/${BUILD_TAG}/refine && \
  ./jenkins/remove_old_builds.sh \
   ${BUILD_NUMBER} \
   "${build_directory}/jenkins-${JOB_NAME}"
EOF
}

