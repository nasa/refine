
module load gcc_6.2.0
module load openmpi_2.1.1_intel_2017
module load intel_2017.2.174

module use --append /u/shared/fun3d/fun3d_users/modulefiles
module load tetgen/1.5.0

module use --append /u/shared/wtjones1/Modules/modulefiles
module load GEOLAB/geolab_64 GEOLAB/AFLR3-16.28.5

module load ESP/119-beta.2020.09.29.1106

module load valgrind_3.13.0

module list

export module_path="/u/shared/fun3d/fun3d_users/modules"

export mpi_path="/usr/local/pkgs-modules/openmpi_2.1.1_intel_2017"

export parmetis_path="${module_path}/ParMETIS/4.0.3-openmpi-2.1.1-intel_2017.2.174"
export zoltan_path="${module_path}/Zoltan/3.82-openmpi-1.10.7-intel_2017.2.174"

export egads_path="${module_path}/ESP/119-beta.2020.09.29.1106/EngSketchPad"
export opencascade_path="${module_path}/ESP/119-beta.2020.09.29.1106/OpenCASCADE-7.3.1"

