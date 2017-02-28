#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

top_dir=$(dirname $PWD)
jenkins_dir=$(dirname $0)
source ${jenkins_dir}/common-jenkins-test-methods

function run_commit_sentinel(){
  echo Running commit sentinel
  ./Ruby/check_fortran.rb `find . -name .git -prune -o -name "*90"`
  ./Ruby/check_makefile_am.rb  `find . -name .git -prune -o -name "*.am"`
  ./Ruby/check_gitignore.rb  `find . -name .git -prune -o -name ".gitignore"`
}

function load_modules() {
  . /usr/local/pkgs/modules/init/bash
  module add gcc_4.9.1_64
  module add make_4.1
  module add Python_2.7.8 # for python build
}

function run_sentinel_quick_start() {
  echo Run sequential quick start
  cd doc/quick_start
      trap "cat step_01.out" EXIT
        env TESTING_PATH=$PATH SRCDIR=. ./make_check.sh
      trap - EXIT
      module add tecplot360-2013R1_64
      ./make_artifacts.sh
      module rm tecplot360-2013R1_64
  cd ../..
}

function compile_fun3d_gfortran() {
  local libcore_location=$1
  local fortran_compiler=$2
  local fc_clags=$3
  local build_dir=$4
  echo "GFortran compilation with real gas phys"
  echo "building with libcore at: $libcore_location"

  if [ $libcore_location == "subpackage" ]; then
    libcore_config_opt=""
    build_dir="_gfortran-subpackage-libcore"
  else
    libcore_config_opt="--with-libcore=${libcore_location}"
    build_dir="_gfortran-external-libcore"
  fi

  mkdir -p $build_dir
  cd $build_dir
    ../configure \
      --prefix=${PWD} \
      ${libcore_config_opt} \
      --enable-complex \
      --enable-hefss \
      --disable-ftune \
      --enable-python \
      FC="${fortran_compiler}" \
      FCFLAGS="${gfortran_fc_flags}" \

    env TMPDIR=${PWD} make -j 8
    env TMPDIR=${PWD} make -j 8 complex
    env TMPDIR=${PWD} make install
    export PATH=${PWD}/bin:$PATH
  cd ..

}

function generate_and_check_manual() {
  echo Extract namelist docs, spell check, and typeset users manual
  cd doc/manual
    ./extract_namelists.rb
    set +x # turn off tracing
    exit_code=0
    for texfile in *.tex ; do
      unknown_words=$( cat $texfile \
                       | aspell -t -p ${PWD}/personal_dictionary.en.pws list \
                       | sort -u )
      if [ -n "$unknown_words" ]; then
        echo "ERROR: ${texfile} has unknown words(s):"
        echo "=>>"
        echo "${unknown_words}"
        echo "<<="
        echo "Please correct spelling and/or add unknown word(s)"
        echo "to private dictionary: doc/manual/personal_dictionary.en.pws"
        echo
        if [ $( echo $texfile | grep '^nml_' ) ] ; then
          echo "Note: This filename has the form 'nml_namelist_name.tex' which"
          echo "was created from Fortran source code where the namelist was"
          echo "defined.  Most live in the libinit, but this should find it:"
          echo "grep -m1 'namelist */namelist_name' \`find . -type f -name \*90\`"
          echo "(where you would replace 'namelist_name' as appropriate)."
          echo
        fi
        exit_code=1
      fi
    done
    set -x # turn tracing back on
    if [ $exit_code -ne 0 ]; then exit $exit_code ; fi
      version=$( grep -oEm1 '[0-9]+\.[0-9]+' ../../configure.ac )
      sed -e "s/M.m/${version}/" manual.tex > FUN3D_Manual-${version}.tex
      export TEXINPUTS="../quick_start:../../:"
      pdflatex FUN3D_Manual-${version}
      bibtex   FUN3D_Manual-${version}
      pdflatex FUN3D_Manual-${version}
      pdflatex FUN3D_Manual-${version}
  cd ../..
}

run_commit_sentinel
load_modules

fortran_compiler="gfortran"
gfortran_fc_flags='-std=f2008 -Werror -Wpedantic -Wall -Wextra -Wsurprising -Wtabs -Wcharacter-truncation -Wunderflow -Wunused-dummy-argument -Wunused-parameter -Walign-commons -Wfunction-elimination -Wzerotrip -fall-intrinsics -Wno-realloc-lhs -Wno-realloc-lhs-all -Wno-array-temporaries -Wno-conversion -Wno-conversion-extra -Wno-compare-reals -Wno-implicit-interface -Wno-implicit-procedure -Wno-aliasing'
config_opts=" --enable-complex --disable-ftune"

libcore_source_dir=`get_libcore_location ${top_dir}`
libcore_install_dir="${libcore_source_dir}/_${fortran_compiler}"
build_libcore \
  "${fortran_compiler}" \
  "${gfortran_fc_flags}" \
  "${libcore_source_dir}" \
  "${libcore_install_dir}" \
  "${config_opts}"

./bootstrap

fortran_compiler="gfortran"
compile_fun3d_gfortran ${libcore_install_dir} ${fortran_compiler} ${gfortran_fc_flags}

compile_fun3d_gfortran "subpackage" ${fortran_compiler} ${gfortran_fc_flags}

run_sentinel_quick_start

generate_and_check_manual


# FIXME: strict compile of PHYSICS_DUMMY
echo GFortran makedist compilation dummy physics
export CONFIG_OPTS="--disable-hefss --disable-ftune \
                    --with-libcore=${libcore_install_dir}"
export FC=gfortran
export FCFLAGS='-O0'
mkdir _gfortran-disable-hefss
cd _gfortran-disable-hefss

    ../configure --prefix=${PWD} $CONFIG_OPTS

    env TMPDIR=${PWD} make -j 8 --output-sync distcheck \
	DISTCHECK_CONFIGURE_FLAGS="$CONFIG_OPTS"

cd ..
unset FC
unset FCFLAGS

echo Build with NAGware compiler
LOG=$PWD/log.nagware
trap "cat $LOG" EXIT
mkdir _nagware
cd _nagware

    module add NAGWare_6.0_64     > $LOG 2>&1

    ../configure \
      --disable-ftune \
      --enable-hefss \
      --without-knife \
      --without-refine \
      FC=nagfor \
      FCFLAGS='-O0 -C=all -f2008' >> $LOG 2>&1

    env TMPDIR=${PWD} make -j 8 --output-sync=target 1>> $LOG 2>&1

    module rm NAGWare_6.0_64
cd ..
trap - EXIT

echo Collect all compiler warnings
grep \
    -e Warning \
    -e Questionable \
    -e Extension \
    -e Obsolescent \
    -e 'Deleted feature used' \
    -e Error \
    -e Fatal \
    -e Panic \
    $LOG | \
    grep \
	-v \
	-e '^make' \
	> all-nagware-warnings || true

echo Search for unexcusable warnings
grep \
    -v \
    -e "Last statement of DO loop body is an unconditional" \
    -e "CONVERT= specifier in OPEN statement" \
    -e "Warning: ignoring incorrect section type for .lbss" \
    all-nagware-warnings \
    > FAILURE-unexcused-nagware-warnings || true
if [ -s FAILURE-unexcused-nagware-warnings ] ; then
    echo FAILURE: Unexcused warnings found
    cat FAILURE-unexcused-nagware-warnings
    exit 1
fi
