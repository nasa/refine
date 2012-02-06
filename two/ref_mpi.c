
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "ref_mpi.h"

REF_INT ref_mpi_n = 1;
REF_INT ref_mpi_id = 0;

REF_STATUS ref_mpi_start( int argc, char *argv[] )
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD,&ref_mpi_n);
  MPI_Comm_rank(MPI_COMM_WORLD,&ref_mpi_id);
#else
  SUPRESS_UNUSED_COMPILER_WARNING(argc);
  SUPRESS_UNUSED_COMPILER_WARNING(argv[0]);
  ref_mpi_n = 1;
  ref_mpi_id = 0;
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_stop( )
{

#ifdef HAVE_MPI
  MPI_Finalize( );
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_bcast( void *data, REF_INT n, REF_TYPE type )
{
#ifdef HAVE_MPI
  MPI_Datatype datatype;

  switch (type)
    {
    case REF_INT_TYPE: datatype = MPI_INT; break;
    case REF_DBL_TYPE: datatype = MPI_DOUBLE; break;
    default: RSS( REF_IMPLEMENT, "data type");
    }

  MPI_Bcast(data, n, datatype, 0, MPI_COMM_WORLD);
#else
  SUPRESS_UNUSED_COMPILER_WARNING(data);
  SUPRESS_UNUSED_COMPILER_WARNING(n);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_send( void *data, REF_INT n, REF_TYPE type, REF_INT dest )
{
#ifdef HAVE_MPI
  MPI_Datatype datatype;
  REF_INT tag;

  switch (type)
    {
    case REF_INT_TYPE: datatype = MPI_INT; break;
    case REF_DBL_TYPE: datatype = MPI_DOUBLE; break;
    default: RSS( REF_IMPLEMENT, "data type");
    }

  tag = ref_mpi_n*dest+ref_mpi_id;

  MPI_Send(data, n, datatype, dest, tag, MPI_COMM_WORLD);
#else
  SUPRESS_UNUSED_COMPILER_WARNING(data);
  SUPRESS_UNUSED_COMPILER_WARNING(n);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
  SUPRESS_UNUSED_COMPILER_WARNING(dest);
  return REF_IMPLEMENT;
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_recv( void *data, REF_INT n, REF_TYPE type, REF_INT source )
{
#ifdef HAVE_MPI
  MPI_Datatype datatype;
  REF_INT tag;
  MPI_Status status;

  switch (type)
    {
    case REF_INT_TYPE: datatype = MPI_INT; break;
    case REF_DBL_TYPE: datatype = MPI_DOUBLE; break;
    default: RSS( REF_IMPLEMENT, "data type");
    }

  tag = ref_mpi_n*ref_mpi_id+source;

  MPI_Recv(data, n, datatype, source, tag, MPI_COMM_WORLD, &status);
#else
  SUPRESS_UNUSED_COMPILER_WARNING(data);
  SUPRESS_UNUSED_COMPILER_WARNING(n);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
  SUPRESS_UNUSED_COMPILER_WARNING(source);
  return REF_IMPLEMENT;
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_alltoall( void *send, void *recv, REF_TYPE type )
{
#ifdef HAVE_MPI
  MPI_Datatype datatype;

  switch (type)
    {
    case REF_INT_TYPE: datatype = MPI_INT; break;
    case REF_DBL_TYPE: datatype = MPI_DOUBLE; break;
    default: RSS( REF_IMPLEMENT, "data type");
    }

  MPI_Alltoall(send, 1, datatype, 
	       recv, 1, datatype, 
	       MPI_COMM_WORLD );
#else
  SUPRESS_UNUSED_COMPILER_WARNING(send);
  SUPRESS_UNUSED_COMPILER_WARNING(recv);
  SUPRESS_UNUSED_COMPILER_WARNING(type);
  return REF_IMPLEMENT;
#endif

  return REF_SUCCESS;
}
