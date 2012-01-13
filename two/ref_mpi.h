
#ifndef REF_MPI_H
#define REF_MPI_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

extern REF_INT ref_mpi_n;
extern REF_INT ref_mpi_id;

#define ref_mpi_master (0 == ref_mpi_id)

REF_STATUS ref_mpi_start( int argc, char *argv[] );
REF_STATUS ref_mpi_stop( );

END_C_DECLORATION

#endif /* REF_MPI_H */
