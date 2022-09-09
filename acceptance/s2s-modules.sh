
esp_module=ESP/121-beta.2022.08.18.1700

module use --append /u/shared/fun3d/fun3d_users/modulefiles
module use --append /u/shared/fun3d/fun3d_users/test_modulefiles

module load --auto FUN3D_INTG
module load ${esp_module}
module load tetgen

export module_path="/u/shared/fun3d/fun3d_users/modules"

export mpi_path="/opt/hpe/hpc/mpt/mpt-2.25"

export parmetis_path="${module_path}/ParMETIS/4.0.3-mpt-2.25-intel_2019.5.281"

export egads_path="${module_path}/${esp_module}/EngSketchPad"
export opencascade_path="${module_path}/${esp_module}/OpenCASCADE"

