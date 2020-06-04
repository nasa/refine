
module load gcc_6.2.0
module load openmpi_1.10.2_intel_2017
module load intel.2017.2.174

module use --append /ump/fldmd/home/casb-shared/fun3d/fun3d_users/modulefiles
module load tetgen/1.5.0

module use --append /usr/local/pkgs-geolab/Modules/modulefiles
module load geolab_64 AFLR3/16.28.5

module load valgrind_3.13.0

module list

export parmetis_path="/ump/fldmd/home/casb-shared/fun3d/fun3d_users/modules/ParMETIS/4.0.3-1.10.2_intel_2017-2017.2.174"
export zoltan_path="/ump/fldmd/home/casb-shared/fun3d/fun3d_users/modules/Zoltan/3.82-1.10.2_intel_2017-2017.2.174"

export egads_path="/ump/fldmd/home/casb-shared/fun3d/fun3d_users/modules/ESP/118/EngSketchPad"
export opencascade_path="/ump/fldmd/home/casb-shared/fun3d/fun3d_users/modules/ESP/114/OpenCASCADE-6.8.1/lin64/gcc"
