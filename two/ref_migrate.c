
#include <stdlib.h>
#include <stdio.h>

#include "ref_migrate.h"

#include "ref_mpi.h"

#ifdef HAVE_ZOLTAN
#include "zoltan.h"
static struct Zoltan_Struct *zz;

static int ref_migrate_local_nodes( void *void_ref_migrate, int *ierr )
{
  REF_MIGRATE ref_migrate = (REF_MIGRATE)void_ref_migrate;
  REF_NODE ref_node = ref_grid_node(ref_migrate_grid(ref_migrate));
  REF_INT node, local_nodes;
  *ierr = 0;
  local_nodes = 0;
  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id == ref_node_part(ref_node,node) ) local_nodes++;
  return local_nodes;
}

static int ref_migrate_geometric_dimensionality( void *void_ref_migrate, 
						 int *ierr )
{
  SUPRESS_UNUSED_COMPILER_WARNING(void_ref_migrate);
  *ierr = 0;
  return 3;
}

#endif

REF_STATUS ref_migrate_create( REF_MIGRATE *ref_migrate_ptr, REF_GRID ref_grid )
{
  REF_MIGRATE ref_migrate;

  (*ref_migrate_ptr) = NULL;
  (*ref_migrate_ptr) = (REF_MIGRATE)malloc( sizeof(REF_MIGRATE_STRUCT) );
  RNS(*ref_migrate_ptr,"malloc ref_migrate NULL");

  ref_migrate = *ref_migrate_ptr;

  ref_migrate_grid(ref_migrate) = ref_grid;

  RSS( ref_edge_create( &(ref_migrate_edge( ref_migrate )), 
			ref_migrate_grid(ref_migrate) ), "create edge" );

#ifdef HAVE_ZOLTAN
#define ref_migrate_zz ((Zoltan_Struct *)ref_migrate->partitioner_data)
  {
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    float ver;
    REIS( ZOLTAN_OK, 
	  Zoltan_Initialize(ref_mpi_argc, ref_mpi_argv, &ver), 
	  "Zoltan is angry");
    zz = Zoltan_Create(MPI_COMM_WORLD);

    /* General parameters */

    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "LB_METHOD", "RCB");

    Zoltan_Set_Num_Obj_Fn(zz, ref_migrate_local_nodes, 
			  (void *)ref_migrate);
    Zoltan_Set_Num_Geom_Fn(zz, ref_migrate_geometric_dimensionality, 
			       (void *)ref_migrate);

    REIS( ZOLTAN_OK, 
	  Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
        &changes,        /* 1 if partitioning was changed, 0 otherwise */ 
        &numGidEntries,  /* Number of integers used for a global ID */
        &numLidEntries,  /* Number of integers used for a local ID */
        &numImport,      /* Number of vertices to be sent to me */
        &importGlobalGids,  /* Global IDs of vertices to be sent to me */
        &importLocalGids,   /* Local IDs of vertices to be sent to me */
        &importProcs,    /* Process rank for source of each incoming vertex */
        &importToPart,   /* New partition for each incoming vertex */
        &numExport,      /* Number of vertices I must send to other processes*/
        &exportGlobalGids,  /* Global IDs of the vertices I must send */
        &exportLocalGids,   /* Local IDs of the vertices I must send */
        &exportProcs,    /* Process to which I send each of the vertices */
        &exportToPart),  /* Partition to which each vertex will belong */
	  "Zoltan is angry");
	  
  }
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_free( REF_MIGRATE ref_migrate )
{
  if ( NULL == (void *)ref_migrate ) return REF_NULL;

#ifdef HAVE_ZOLTAN
  Zoltan_Destroy( &zz );
#endif

  RSS( ref_edge_free( ref_migrate_edge( ref_migrate ) ), "free edge" );

  ref_cond_free( ref_migrate );

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_inspect( REF_MIGRATE ref_migrate )
{
  RSS(ref_grid_inspect( ref_migrate_grid(ref_migrate)) , "inspect grid");

  return REF_SUCCESS;
}

