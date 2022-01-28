

module use --append /u/shared/fun3d/fun3d_users/modulefiles
module use --append /u/shared/wtjones1/Modules/modulefiles

export compile_modules="gcc_6.2.0 intel_2017.2.174"
export compiler_rpath="/usr/local/pkgs-modules/gcc_6.2.0/lib64"

export run_modules="openmpi_2.1.1_intel_2017 tetgen ESP/120 GEOLAB/geolab_64 GEOLAB/AFLR3-16.28.5 valgrind_3.13.0"

module load ${compile_modules} ${run_modules}

export module_path="/u/shared/fun3d/fun3d_users/modules"

export mpi_path="/usr/local/pkgs-modules/openmpi_2.1.1_intel_2017"

export parmetis_path="${module_path}/ParMETIS/4.0.3-openmpi-2.1.1-intel_2017.2.174"
export zoltan_path="${module_path}/Zoltan/3.82-openmpi-1.10.7-intel_2017.2.174"

export egads_path="${module_path}/ESP/120/EngSketchPad"
export opencascade_path="${module_path}/ESP/120/OpenCASCADE"
