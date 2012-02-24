
#include <stdlib.h>
#include <stdio.h>

#include "ref_migrate.h"

#include "ref_mpi.h"

#ifdef HAVE_ZOLTAN
#include "zoltan.h"
static struct Zoltan_Struct *zz;

static int ref_migrate_local_n( void *void_ref_migrate, int *ierr )
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

static void ref_migrate_local_ids( void *void_ref_migrate, 
				   int global_dim, int local_dim,
				   ZOLTAN_ID_PTR global, ZOLTAN_ID_PTR local,
				   int wgt_dim, float *obj_wgts, int *ierr )
{
  REF_MIGRATE ref_migrate = (REF_MIGRATE)void_ref_migrate;
  REF_NODE ref_node = ref_grid_node(ref_migrate_grid(ref_migrate));
  REF_INT node, nnode;
  if ( 1 != global_dim || 1 != local_dim || 0 !=  wgt_dim )
    {
      printf("%s: %d: %s: %s\n",__FILE__,__LINE__,__func__,"bad sizes");
      *ierr = ZOLTAN_FATAL;
      return;
    }
  SUPRESS_UNUSED_COMPILER_WARNING(obj_wgts);
  *ierr = 0;
  nnode = 0;
  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id == ref_node_part(ref_node,node) ) 
      {
	local[node] = node;
	global[node] = ref_node_global(ref_node,node);
	nnode++;
      }
}

static int ref_migrate_geometric_dimensionality( void *void_ref_migrate, 
						 int *ierr )
{
  SUPRESS_UNUSED_COMPILER_WARNING(void_ref_migrate);
  *ierr = 0;
  return 3;
}

static void ref_migrate_xyz( void *void_ref_migrate, 
			     int global_dim, int local_dim, int nnode,
			     ZOLTAN_ID_PTR global, ZOLTAN_ID_PTR local,
			     int xyz_dim, double *xyz, int *ierr )
{
  REF_MIGRATE ref_migrate = (REF_MIGRATE)void_ref_migrate;
  REF_NODE ref_node = ref_grid_node(ref_migrate_grid(ref_migrate));
  REF_INT node;
  if ( 1 != global_dim || 1 != local_dim || 3 !=  xyz_dim )
    {
      printf("%s: %d: %s: %s\n",__FILE__,__LINE__,__func__,"bad sizes");
      *ierr = ZOLTAN_FATAL;
      return;
    }
  SUPRESS_UNUSED_COMPILER_WARNING(global);
  *ierr = 0;
  for (node=0;node<nnode;node++)
    {
      xyz[0+3*node] = ref_node_xyz(ref_node,0,local[node]);
      xyz[1+3*node] = ref_node_xyz(ref_node,1,local[node]);
      xyz[2+3*node] = ref_node_xyz(ref_node,2,local[node]);
    }

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

#ifdef HAVE_ZOLTAN
#define ref_migrate_zz ((Zoltan_Struct *)ref_migrate->partitioner_data)
  {
    int partitions_have_changed;
    int global_id_dimension, local_id_dimension;

    int import_n;
    ZOLTAN_ID_PTR import_global, import_local;
    int *import_proc, *import_part;

    int export_n;
    ZOLTAN_ID_PTR export_global, export_local;
    int *export_proc, *export_part;

    float ver;
    REIS( ZOLTAN_OK, 
	  Zoltan_Initialize(ref_mpi_argc, ref_mpi_argv, &ver), 
	  "Zoltan is angry");
    zz = Zoltan_Create(MPI_COMM_WORLD);

    /* General parameters */

    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTS");
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "LB_METHOD", "RCB");

    Zoltan_Set_Num_Obj_Fn(zz, ref_migrate_local_n, 
			  (void *)ref_migrate);
    Zoltan_Set_Obj_List_Fn(zz, ref_migrate_local_ids, 
			   (void *)ref_migrate);
    Zoltan_Set_Num_Geom_Fn(zz, ref_migrate_geometric_dimensionality, 
			   (void *)ref_migrate);
    Zoltan_Set_Geom_Multi_Fn(zz, ref_migrate_xyz, 
			     (void *)ref_migrate);

    REIS( ZOLTAN_OK, 
	  Zoltan_LB_Partition(zz,
			      &partitions_have_changed,
			      &global_id_dimension,
			      &local_id_dimension,
			      &import_n,
			      &import_global,
			      &import_local,
			      &import_proc,
			      &import_part,
			      &export_n,
			      &export_global,
			      &export_local,
			      &export_proc,
			      &export_part),
	  "Zoltan is angry");

    REIS( ZOLTAN_OK,
	  Zoltan_LB_Free_Part(&import_local, &import_global,
			      &import_proc, &import_proc ),
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

  ref_cond_free( ref_migrate );

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_inspect( REF_MIGRATE ref_migrate )
{
  RSS(ref_grid_inspect( ref_migrate_grid(ref_migrate)) , "inspect grid");

  return REF_SUCCESS;
}

