
#ifndef REF_MPI_H
#define REF_MPI_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

extern REF_INT ref_mpi_n;
extern REF_INT ref_mpi_id;

#define ref_mpi_master (0 == ref_mpi_id)

typedef int REF_TYPE;
#define REF_INT_TYPE (1)
#define REF_DBL_TYPE (2)

REF_STATUS ref_mpi_start( int argc, char *argv[] );
REF_STATUS ref_mpi_stop( );

REF_STATUS ref_mpi_bcast( void *data, REF_INT n, REF_TYPE type );

REF_STATUS ref_mpi_send( void *data, REF_INT n, REF_TYPE type, REF_INT dest );
REF_STATUS ref_mpi_recv( void *data, REF_INT n, REF_TYPE type, REF_INT source );

REF_STATUS ref_mpi_alltoall( void *send, void *recv, REF_TYPE type );
REF_STATUS ref_mpi_alltoallv( void *send, REF_INT *send_size, 
			      void *recv, REF_INT *recv_size, 
			      REF_INT n, REF_TYPE type );

END_C_DECLORATION

#endif /* REF_MPI_H */
