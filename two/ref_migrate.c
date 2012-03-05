
#include <stdlib.h>
#include <stdio.h>

#include "ref_migrate.h"

#include "ref_mpi.h"
#include "ref_malloc.h"

#include "ref_part.h"

#include "ref_sort.h"
#include "ref_node.h"

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

static REF_STATUS ref_migrate_shufflin_node( REF_NODE ref_node )
{
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT *a_global, *b_global;
  REF_INT part, node;
  REF_INT *a_next;
  REF_DBL *a_xyz, *b_xyz;
  REF_INT local;

  ref_malloc_init( a_size, ref_mpi_n, REF_INT, 0 );
  ref_malloc_init( b_size, ref_mpi_n, REF_INT, 0 );

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      a_size[ref_node_part(ref_node,node)]++;

  RSS( ref_mpi_alltoall( a_size, b_size, REF_INT_TYPE ), "alltoall sizes");

  a_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    a_total += a_size[part];
  ref_malloc( a_global, a_total, REF_INT );
  ref_malloc( a_xyz, 3*a_total, REF_DBL );

  b_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    b_total += b_size[part];
  ref_malloc( b_global, b_total, REF_INT );
  ref_malloc( b_xyz, 3*b_total, REF_DBL );

  ref_malloc( a_next, ref_mpi_n, REF_INT );
  a_next[0] = 0;
  for ( part = 1; part<ref_mpi_n ; part++ )
    a_next[part] = a_next[part-1]+a_size[part-1];

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      {
	part = ref_node_part(ref_node,node);
	a_global[a_next[part]] = ref_node_global(ref_node,node);
	a_xyz[0+3*a_next[part]] = ref_node_xyz(ref_node,0,node);
	a_xyz[1+3*a_next[part]] = ref_node_xyz(ref_node,1,node);
	a_xyz[2+3*a_next[part]] = ref_node_xyz(ref_node,2,node);
	a_next[part]++;
      }

  RSS( ref_mpi_alltoallv( a_global, a_size, b_global, b_size, 
			  1, REF_INT_TYPE ), 
       "alltoallv global");

  RSS( ref_mpi_alltoallv( a_xyz, a_size, b_xyz, b_size, 
			  3, REF_DBL_TYPE ), 
       "alltoallv global");

  RSS( ref_node_add_many( ref_node, b_total, b_global ), "add many" );

  for ( node=0; node < b_total; node++ )
    {
      RSS( ref_node_local( ref_node, b_global[node], &local ), "local" );
      ref_node_xyz(ref_node,0,local) = b_xyz[0+3*node];
      ref_node_xyz(ref_node,1,local) = b_xyz[1+3*node];
      ref_node_xyz(ref_node,2,local) = b_xyz[2+3*node];
      ref_node_part(ref_node,local) = ref_mpi_id;
    }

  free(a_next);
  free(b_xyz);
  free(b_global);
  free(a_xyz);
  free(a_global);
  free(b_size);
  free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_shufflin_cell( REF_NODE ref_node, 
				      REF_CELL ref_cell )
{
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT all_parts[REF_CELL_MAX_SIZE_PER];
  REF_INT nunique;
  REF_INT unique_parts[REF_CELL_MAX_SIZE_PER];
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT part, node, cell, i;
  REF_INT *a_next;
  REF_INT *a_c2n, *b_c2n;
  REF_INT *a_parts, *b_parts;
  REF_BOOL need_to_keep;

  if ( 1 == ref_mpi_n ) return REF_SUCCESS;

  ref_malloc_init( a_size, ref_mpi_n, REF_INT, 0 );
  ref_malloc_init( b_size, ref_mpi_n, REF_INT, 0 );

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      for ( node=0; node < ref_cell_node_per(ref_cell); node++ )
	all_parts[node] = ref_node_part(ref_node,nodes[node]);
      RSS( ref_sort_unique( ref_cell_node_per(ref_cell), all_parts,
			    &nunique, unique_parts ), "unique");
      for ( node=0; node < nunique; node++ )
	{
	  part = unique_parts[node];
	  if ( ref_mpi_id != part ) a_size[part]++;
	}
    }

  RSS( ref_mpi_alltoall( a_size, b_size, REF_INT_TYPE ), "alltoall sizes");

  a_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    a_total += a_size[part];
  ref_malloc( a_c2n,   ref_cell_size_per(ref_cell)*a_total, REF_INT );
  ref_malloc( a_parts, ref_cell_size_per(ref_cell)*a_total, REF_INT );

  b_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    b_total += b_size[part];
  ref_malloc( b_c2n,   ref_cell_size_per(ref_cell)*b_total, REF_INT );
  ref_malloc( b_parts, ref_cell_size_per(ref_cell)*b_total, REF_INT );

  ref_malloc( a_next, ref_mpi_n, REF_INT );
  a_next[0] = 0;
  for ( part = 1; part<ref_mpi_n ; part++ )
    a_next[part] = a_next[part-1]+a_size[part-1];

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      for ( node=0; node < ref_cell_node_per(ref_cell); node++ )
	all_parts[node] = ref_node_part(ref_node,nodes[node]);
      RSS( ref_sort_unique( ref_cell_node_per(ref_cell), all_parts,
			    &nunique, unique_parts ), "unique");
      for ( node=0; node < nunique; node++ )
	{
	  part = unique_parts[node];
	  if ( ref_mpi_id != part ) 
	    {
	      for (i=0;i<ref_cell_node_per(ref_cell);i++)
		{
		  a_c2n[i+ref_cell_size_per(ref_cell)*a_next[part]] = 
		    ref_node_global(ref_node,nodes[i]);
		  a_parts[i+ref_cell_size_per(ref_cell)*a_next[part]] =
		    ref_node_part(ref_node,nodes[i]);
		}
	      if ( ref_cell_last_node_is_an_id(ref_cell) )
		a_c2n[ref_cell_node_per(ref_cell) + 
		      ref_cell_size_per(ref_cell)*a_next[part]] =
		  nodes[ref_cell_node_per(ref_cell)];
	      a_next[part]++;
	    }
	}
    }

  RSS( ref_mpi_alltoallv( a_c2n, a_size, b_c2n, b_size, 
			  ref_cell_size_per(ref_cell), REF_INT_TYPE ), 
       "alltoallv c2n");
  RSS( ref_mpi_alltoallv( a_parts, a_size, b_parts, b_size, 
			  ref_cell_size_per(ref_cell), REF_INT_TYPE ), 
       "alltoallv parts");

  RSS( ref_cell_add_many_global( ref_cell, ref_node,
				 b_total, 
				 b_c2n, b_parts, ref_mpi_id ), "many glob");

  free(a_next);
  free(b_c2n);
  free(a_c2n);
  free(b_size);
  free(a_size);

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      need_to_keep = REF_FALSE;
      for ( node=0; node < ref_cell_node_per(ref_cell); node++ )
	{
	  need_to_keep = ( need_to_keep || 
			   ref_mpi_id == ref_node_part(ref_node,nodes[node]) );
	}
      if ( ! need_to_keep )
	RSS( ref_cell_remove( ref_cell, cell), "remove" );
    }

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_shufflin( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node( ref_grid );
  REF_CELL ref_cell;
  REF_INT group, node;
  REF_BOOL need_to_keep;

  if ( 1 == ref_mpi_n ) return REF_SUCCESS;

  ref_mpi_stopwatch_start();
  RSS( ref_migrate_shufflin_node( ref_node ), "send out nodes" );
  ref_mpi_stopwatch_stop("node");

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    {
      RSS( ref_migrate_shufflin_cell( ref_node, ref_cell ), "cell" );
      ref_mpi_stopwatch_stop("cell");
    }

  RSS( ref_migrate_shufflin_cell( ref_node, ref_grid_tri(ref_grid) ), "tri");
  ref_mpi_stopwatch_stop("tri");
  RSS( ref_migrate_shufflin_cell( ref_node, ref_grid_qua(ref_grid) ), "qua");
  ref_mpi_stopwatch_stop("qua");

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      {
	need_to_keep = REF_FALSE;
	each_ref_grid_ref_cell( ref_grid, group, ref_cell )
	  need_to_keep = ( need_to_keep ||
			   !ref_adj_empty( ref_cell_adj(ref_cell), node ) );
	if ( !need_to_keep )
	  RSS( ref_node_remove_without_global( ref_node, node), "remove" );
      }
  RSS( ref_node_rebuild_sorted_global( ref_node ), "rebuild" );

  ref_mpi_stopwatch_stop("remove");

  RSS( ref_part_ghost_xyz( ref_grid ), "ghost xyz");

  ref_mpi_stopwatch_stop("ghost");

  return REF_SUCCESS;
}

