
#include <stdlib.h>
#include <stdio.h>

#include "ref_migrate.h"

#include "ref_mpi.h"
#include "ref_malloc.h"

#include "ref_part.h"

#include "ref_export.h"

#ifdef HAVE_ZOLTAN
#include "zoltan.h"

static int ref_migrate_local_n( void *void_ref_grid, int *ierr )
{
  REF_NODE ref_node = ref_grid_node((REF_GRID)void_ref_grid);
  REF_INT node, local_nodes;
  *ierr = 0;
  local_nodes = 0;
  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id == ref_node_part(ref_node,node) ) local_nodes++;
  return local_nodes;
}

static void ref_migrate_local_ids( void *void_ref_grid, 
				   int global_dim, int local_dim,
				   ZOLTAN_ID_PTR global, ZOLTAN_ID_PTR local,
				   int wgt_dim, float *obj_wgts, int *ierr )
{
  REF_NODE ref_node = ref_grid_node((REF_GRID)void_ref_grid);
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

static int ref_migrate_geometric_dimensionality( void *void_ref_grid, 
						 int *ierr )
{
  SUPRESS_UNUSED_COMPILER_WARNING(void_ref_grid);
  *ierr = 0;
  return 3;
}

static void ref_migrate_xyz( void *void_ref_grid, 
			     int global_dim, int local_dim, int nnode,
			     ZOLTAN_ID_PTR global, ZOLTAN_ID_PTR local,
			     int xyz_dim, double *xyz, int *ierr )
{
  REF_NODE ref_node = ref_grid_node((REF_GRID)void_ref_grid);
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

REF_STATUS ref_migrate_new_part( REF_GRID ref_grid )
{

#ifdef HAVE_ZOLTAN
  {
    REF_NODE ref_node = ref_grid_node( ref_grid );
    int partitions_have_changed;
    int global_id_dimension, local_id_dimension;

    int import_n;
    ZOLTAN_ID_PTR import_global, import_local;
    int *import_proc, *import_part;

    int export_n;
    ZOLTAN_ID_PTR export_global, export_local;
    int *export_proc, *export_part;

    float ver;

    REF_INT node;

    REF_INT *part;

    struct Zoltan_Struct *zz;

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
			  (void *)ref_grid);
    Zoltan_Set_Obj_List_Fn(zz, ref_migrate_local_ids, 
			   (void *)ref_grid);
    Zoltan_Set_Num_Geom_Fn(zz, ref_migrate_geometric_dimensionality, 
			   (void *)ref_grid);
    Zoltan_Set_Geom_Multi_Fn(zz, ref_migrate_xyz, 
			     (void *)ref_grid);

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

    ref_malloc_init( part, ref_node_max(ref_node), REF_INT, REF_EMPTY );

    for(node=0; node<export_n; node++)
      part[export_local[node]] = export_part[node];

    RSS( ref_part_ghost_int( ref_grid, part ), "ghost part");

    for(node=0; node<ref_node_max(ref_node); node++)
      ref_node_part(ref_node, node) = part[node];

    ref_free( part );

    REIS( ZOLTAN_OK,
	  Zoltan_LB_Free_Part(&import_local, &import_global,
			      &import_proc, &import_proc ),
	  "Zoltan is angry");

    REIS( ZOLTAN_OK,
	  Zoltan_LB_Free_Part(&export_local, &export_global,
			      &export_proc, &export_proc ),
	  "Zoltan is angry");

  }
#else
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_part_viz( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node( ref_grid );
  char viz_file[256];

  sprintf(viz_file, "ref_migrate_n%d_p%d.tec", ref_mpi_n, ref_mpi_id);

  RSS(ref_export_tec_int( ref_grid, ref_node->part,
			  viz_file ) , "viz parts as scalar");

  return REF_SUCCESS;
}

