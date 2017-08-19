
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_EGADS
#include "egads.h"
#endif

#include "ref_geom.h"

#include "ref_export.h"
#include "ref_gather.h"
#include "ref_grid.h"
#include "ref_node.h"
#include "ref_cell.h"

#include "ref_dict.h"

#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_mpi.h"

REF_STATUS ref_geom_create( REF_GEOM *ref_geom_ptr )
{
  REF_GEOM ref_geom;
  REF_INT geom;
  ( *ref_geom_ptr ) = NULL;

  ref_malloc( *ref_geom_ptr, 1, REF_GEOM_STRUCT );

  ref_geom = ( *ref_geom_ptr );

  ref_geom_n(ref_geom) = 0;
  ref_geom_max(ref_geom) = 10;

  ref_malloc( ref_geom->descr, 3*ref_geom_max(ref_geom), REF_INT);
  ref_malloc( ref_geom->param, 2*ref_geom_max(ref_geom), REF_DBL);
  ref_geom->uv_area_sign = NULL;
  
  for ( geom = 0; geom < ref_geom_max(ref_geom); geom++ )
    {
      ref_geom_type(ref_geom,geom) = REF_EMPTY;
      ref_geom_id(ref_geom,geom) = geom+1;
    }
  ref_geom_id(ref_geom,ref_geom_max(ref_geom)-1) = REF_EMPTY;
  ref_geom_blank(ref_geom) = 0;
  
  RSS( ref_adj_create( &( ref_geom->ref_adj ) ), "create ref_adj for ref_geom" );
  ref_geom->nnode = REF_EMPTY;
  ref_geom->nedge = REF_EMPTY;
  ref_geom->nface = REF_EMPTY;
  ref_geom->context = NULL;
#ifdef HAVE_EGADS
  {
    ego context;
    REIS( EGADS_SUCCESS, EG_open(&context), "EG open");
    /* Success returns the old output level. (0-silent to 3-debug) */
    RAS( EG_setOutLevel(context, 0) >= 0, "make silent");
    ref_geom->context = (void *)context;
  }
#endif
  ref_geom->solid = NULL;
  ref_geom->faces = NULL;
  ref_geom->edges = NULL;
  ref_geom->nodes = NULL;

  ref_geom->cad_data_size = 0;
  ref_geom->cad_data = (REF_BYTE *)NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_free( REF_GEOM ref_geom )
{
  if ( NULL == (void *)ref_geom )
    return REF_NULL;
  ref_free( ref_geom->cad_data );
#ifdef HAVE_EGADS
  if ( NULL != ref_geom->faces)
    EG_free((ego *)(ref_geom->faces));
  if ( NULL != ref_geom->edges)
    EG_free((ego *)(ref_geom->edges));
  if ( NULL != ref_geom->nodes)
    EG_free((ego *)(ref_geom->nodes));
  /* solid is not freeable */
  if ( NULL != ref_geom->context)
    REIS( EGADS_SUCCESS, EG_close((ego)(ref_geom->context)), "EG close");
#endif
  RSS( ref_adj_free( ref_geom->ref_adj ), "adj free" );
  ref_free( ref_geom->uv_area_sign );
  ref_free( ref_geom->param );
  ref_free( ref_geom->descr );
  ref_free( ref_geom );
  return REF_SUCCESS;
}

REF_STATUS ref_geom_deep_copy( REF_GEOM *ref_geom_ptr, REF_GEOM original )
{
  REF_GEOM ref_geom;
  REF_INT geom, i;
  ( *ref_geom_ptr ) = NULL;

  ref_malloc( *ref_geom_ptr, 1, REF_GEOM_STRUCT );

  ref_geom = ( *ref_geom_ptr );

  ref_geom_n(ref_geom) = ref_geom_n(original);
  ref_geom_max(ref_geom) = ref_geom_max(original);

  ref_malloc( ref_geom->descr, 3*ref_geom_max(ref_geom), REF_INT);
  ref_malloc( ref_geom->param, 2*ref_geom_max(ref_geom), REF_DBL);
  ref_geom->uv_area_sign = NULL;
 
  for ( geom = 0; geom < ref_geom_max(ref_geom); geom++ )
    for ( i = 0; i < 3; i++ )
      ref_geom_descr(ref_geom,i,geom) = ref_geom_descr(original,i,geom);
  ref_geom_blank(ref_geom) = ref_geom_blank(original);
  for ( geom = 0; geom < ref_geom_max(ref_geom); geom++ )
    for ( i = 0; i < 2; i++ )
      ref_geom_param(ref_geom,i,geom) = ref_geom_param(original,i,geom);

  RSS( ref_adj_deep_copy( &( ref_geom->ref_adj ), original->ref_adj ),
       "deep copy ref_adj for ref_geom" );
  
  ref_geom->nnode = REF_EMPTY;
  ref_geom->nedge = REF_EMPTY;
  ref_geom->nface = REF_EMPTY;
  ref_geom->context = NULL;
  ref_geom->solid = NULL;
  ref_geom->faces = NULL;
  ref_geom->edges = NULL;
  ref_geom->nodes = NULL;

  ref_geom->cad_data_size = 0;
  ref_geom->cad_data = (REF_BYTE *)NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_save( REF_GRID ref_grid, const char *filename )
{
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_edg(ref_grid);
  FILE *file;
  REF_INT geom, global;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

  RSS( ref_node_synchronize_globals( ref_node ), "sync" );
  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );
  fprintf(file,"%d\n", ref_cell_n(ref_cell) );
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      nodes[0] = ref_node_global(ref_node,nodes[0]);
      nodes[1] = ref_node_global(ref_node,nodes[1]);
      fprintf(file,"%d %d %d\n", nodes[0], nodes[1], nodes[2]);
    }
  fprintf(file,"%d\n", ref_geom_n(ref_geom) );
  each_ref_geom( ref_geom, geom )
    {
      global = ref_node_global(ref_node,ref_geom_node(ref_geom,geom));
      fprintf(file,"%d %d %d %.18e %.18e\n",
	      ref_geom_type(ref_geom,geom),
	      ref_geom_id(ref_geom,geom),
	      global,
	      ref_geom_param(ref_geom,0,geom),
	      ref_geom_param(ref_geom,1,geom));
    }
  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_load( REF_GRID ref_grid, const char *filename )
{
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  FILE *file;
  REF_INT edge, nedge;
  REF_INT new_cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT geom, ngeom;
  REF_INT global, node, gtype, id;
  REF_DBL param[2];
  RSS( ref_node_synchronize_globals( ref_node ), "sync" );
  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );
  REIS( 1, fscanf(file, "%d", &nedge ), "nedge");
  for (edge=0;edge<nedge;edge++)
    {
      REIS( 3, fscanf(file, "%d %d %d",
		      &(nodes[0]), &(nodes[1]), &(nodes[2]) ), "edge nodes");
      RSS( ref_node_local( ref_node, nodes[0], &node ), "g2l" );
      nodes[0] = node;
      RSS( ref_node_local( ref_node, nodes[1], &node ), "g2l" );
      nodes[1] = node;
      RSS( ref_cell_add(ref_grid_edg(ref_grid), nodes, &new_cell ), "new edge");
    }
  REIS( 1, fscanf(file, "%d", &ngeom ), "ngeom");
  for (geom=0;geom<ngeom;geom++)
    {
      REIS( 5, fscanf( file, "%d %d %d %lf %lf\n",
		       &gtype, &id, &global, &(param[0]), &(param[1]) ), "gr" );
      RSS( ref_node_local( ref_node, global, &node ), "g2l" );
      RSS( ref_geom_add( ref_geom, node, gtype, id, param), "add geom");
    }
  if (!feof(file)) 
    {
      printf("expected end of %s file after %d edge %d geom\n",
	     filename, nedge, ngeom );
      fclose(file);
      return REF_FAILURE;
    }
  fclose(file);
  RSS( ref_geom_uv_area_report( ref_grid ), "report");
  return REF_SUCCESS;
}

REF_STATUS ref_geom_recon( REF_GRID ref_grid )
{
#ifdef HAVE_EGADS
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  ego ref, *pchldrn, object;
  int oclass, mtype, nchild, *psens;
  double xyz[9], trange[2];
  REF_INT node, id, best_node;
  REF_DBL best_dist, dist;
  REF_INT *tessnodes;
  REF_INT degree, max_node = 50, *node_list;
  ref_malloc(node_list,max_node,REF_INT);
  printf("searching for %d topo nodes\n",ref_geom->nnode);
  ref_malloc(tessnodes,ref_geom->nnode,REF_INT);
  for ( id = 1 ; id <= ref_geom->nnode ; id++ )
    {
      object = ((ego *)(ref_geom->nodes))[id - 1]; 
      REIS( EGADS_SUCCESS,
	    EG_getTopology(object, &ref, &oclass, &mtype, xyz,
			   &nchild, &pchldrn, &psens), "EG topo node");
      best_node = REF_EMPTY;
      best_dist = 1.0e20;
      each_ref_node_valid_node( ref_node, node )
	{
	  dist = sqrt( pow(xyz[0]-ref_node_xyz(ref_node,0,node),2) +
		       pow(xyz[1]-ref_node_xyz(ref_node,1,node),2) +
		       pow(xyz[2]-ref_node_xyz(ref_node,2,node),2) );
	  if ( dist < best_dist )
	    {
	      best_node = node;
	      best_dist = dist;
	    }
	}
      printf(" topo node id %3d node %6d dist %.4e\n",
	     id,best_node,best_dist);
      tessnodes[id-1]=best_node;
      RSS( ref_geom_add( ref_geom, best_node, REF_GEOM_NODE, id, NULL), "node");
    }
  for ( id = 1 ; id <= ref_geom->nedge ; id++ )
    {
      object = ((ego *)(ref_geom->edges))[id - 1]; 
      REIS( EGADS_SUCCESS,
	    EG_getTopology(object, &ref, &oclass, &mtype, trange,
			   &nchild, &pchldrn, &psens), "EG topo node");
      if (mtype == DEGENERATE)
	{
	  printf(" topo edge id %3d degen\n",id);
	}
      else
	{
	  int toponode0, toponode1;
	  REF_INT node0, node1;
	  double t;
	  double param[2], closest[3];
	  REF_INT i, next_node, current_node;
	  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
	  REF_INT cell, geom;
	  REIS( TWONODE, mtype, "ONENODE edge not implemented");
	  REIS( 2, nchild, "expect to topo node for edge");
	  toponode0 = EG_indexBodyTopo(ref_geom->solid, pchldrn[0]);
	  toponode1 = EG_indexBodyTopo(ref_geom->solid, pchldrn[1]);
	  node0 = tessnodes[toponode0-1];
	  node1 = tessnodes[toponode1-1];
	  printf(" topo edge id %3d\n",id);
	  REIS( EGADS_SUCCESS,
		EG_evaluate(object, &(trange[0]), xyz ), "EG eval");
	  node=node0;
	  dist = sqrt( pow(xyz[0]-ref_node_xyz(ref_node,0,node),2) +
		       pow(xyz[1]-ref_node_xyz(ref_node,1,node),2) +
		       pow(xyz[2]-ref_node_xyz(ref_node,2,node),2) );
	  printf("  node0 id %2d index %3d t %f dist %e\n",
		 toponode0,node0,trange[0],dist);
	  REIS( EGADS_SUCCESS,
		EG_evaluate(object, &(trange[1]), xyz ), "EG eval");
	  node=node1;
	  dist = sqrt( pow(xyz[0]-ref_node_xyz(ref_node,0,node),2) +
		       pow(xyz[1]-ref_node_xyz(ref_node,1,node),2) +
		       pow(xyz[2]-ref_node_xyz(ref_node,2,node),2) );
	  printf("  node1 id %2d index %3d t %f dist %e\n",
		 toponode1,node1,trange[1],dist);
	  current_node = node0;
	  t = trange[0];
	  while ( current_node != node1 )
	    {
	      RSS(ref_cell_node_list_around( ref_grid_tri(ref_grid),
					     current_node, max_node,
					     &degree, node_list),
		  "next node" );
	      best_node = REF_EMPTY;
	      best_dist = 1.0e20;
	      for (i=0;i<degree;i++)
		{
		  next_node = node_list[i];
		  if ( REF_SUCCESS == ref_geom_find( ref_geom, next_node,
						     REF_GEOM_EDGE, id, &geom ))
		    continue; /* this canidate is already part of the edge */
		  REIS( EGADS_SUCCESS,
			EG_invEvaluate(object,
				       ref_node_xyz_ptr(ref_node,next_node),
				       param, closest), "EG eval");
		  dist = sqrt( pow(closest[0]-ref_node_xyz(ref_node,0,next_node),2) +
			       pow(closest[1]-ref_node_xyz(ref_node,1,next_node),2) +
			       pow(closest[2]-ref_node_xyz(ref_node,2,next_node),2) );
		  if ( dist < best_dist )
		    {
		      best_node = next_node;
		      best_dist = dist; 
		    }
		}
	      REIS( EGADS_SUCCESS,
		    EG_invEvaluate(object,
				   ref_node_xyz_ptr(ref_node,best_node),
				   param, closest), "EG eval");
	      printf("   best_node %5d t %f best_dist %e\n",
		     best_node,param[0],best_dist);
	      param[1] = t;
	      nodes[0] = best_node;
	      nodes[1] = current_node;
	      nodes[2] = id;
	      RSS(ref_cell_add(ref_grid_edg(ref_grid), nodes, &cell ), "add e");
	      RSS(ref_geom_add(ref_geom, nodes[0], REF_GEOM_EDGE,
			       id, &(param[0]) ), "add geom edge 0");
	      RSS(ref_geom_add(ref_geom, nodes[1], REF_GEOM_EDGE,
			       id, &(param[1]) ), "add geom edge 1");
	      current_node = best_node;
	      t = param[0];
	    }
	}
    }
  ref_free(tessnodes);
  ref_free(node_list);
  return REF_SUCCESS;
#else
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  if (REF_EMPTY == ref_geom->nnode)
    printf("No EGADS loaded\n");
  printf("No EGADS linked for %s\n", __func__);
  return REF_IMPLEMENT;
#endif

}

REF_STATUS ref_geom_uv_area( REF_GEOM ref_geom, REF_INT *nodes,
			     REF_DBL *uv_area )
{
  REF_INT id = nodes[3];
  REF_DBL uv0[2], uv1[2], uv2[2];
  RSS( ref_geom_tuv( ref_geom, nodes[0], REF_GEOM_FACE, id, uv0 ), "uv0" );
  RSS( ref_geom_tuv( ref_geom, nodes[1], REF_GEOM_FACE, id, uv1 ), "uv1" );
  RSS( ref_geom_tuv( ref_geom, nodes[2], REF_GEOM_FACE, id, uv2 ), "uv2" );
  *uv_area = 0.5 * ( -uv1[0]*uv0[1] + uv2[0]*uv0[1] + uv0[0]*uv1[1]
		     -uv2[0]*uv1[1] - uv0[0]*uv2[1] + uv1[0]*uv2[1] );
  return REF_SUCCESS;
}

REF_STATUS ref_geom_uv_area_sign( REF_GRID ref_grid, REF_INT id,
				  REF_DBL *sign )
{
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  if ( NULL == ((ref_geom)->uv_area_sign) )
    {
      REF_CELL ref_cell = ref_grid_tri(ref_grid);
      REF_INT face;
      REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
      REF_DBL uv_area;
      if (REF_EMPTY == ref_geom->nface)
	RSS( ref_geom_infer_nedge_nface( ref_grid ), "infer counts" );
      ref_malloc_init( ref_geom->uv_area_sign, ref_geom->nface, REF_DBL, 0.0);
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  face = nodes[3];
	  if ( face < 1 || ref_geom->nface < face ) continue;
	  RSS( ref_geom_uv_area( ref_geom, nodes, &uv_area), "uv area");
	  ((ref_geom)->uv_area_sign)[face-1] += uv_area;
	}
      for (face=0;face<ref_geom->nface ;face++)
	{
	  if ( ((ref_geom)->uv_area_sign)[face] < 0.0 )
	    {
	      ((ref_geom)->uv_area_sign)[face] = -1.0;
	    }
	  else
	    {
	      ((ref_geom)->uv_area_sign)[face] = 1.0;
	    }
	}
    }

  if ( id < 1 || id > ref_geom->nface ) return REF_INVALID;
  *sign = ((ref_geom)->uv_area_sign)[id-1];
  
  return REF_SUCCESS;
}

REF_STATUS ref_geom_uv_area_report( REF_GRID ref_grid )
{
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT geom, id, min_id, max_id;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL no_cell;
  REF_DBL uv_area, total_uv_area, min_uv_area, max_uv_area, sign_uv_area;
  REF_INT n_neg, n_pos;
  
  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_face( ref_geom, geom )
    {
      min_id = MIN( min_id, ref_geom_id(ref_geom,geom) );
      max_id = MAX( max_id, ref_geom_id(ref_geom,geom) );
    }

  for ( id = min_id ; id <= max_id ; id++ )
    {
      no_cell = REF_TRUE;
      total_uv_area = 0.0;
      min_uv_area = 0.0;
      max_uv_area = 0.0;
      n_neg = 0;
      n_pos = 0;
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	if ( id == nodes[3] )
	  {
	    RSS( ref_geom_uv_area( ref_geom, nodes, &uv_area), "uv area");
	    total_uv_area += uv_area;
	    if (no_cell)
	      {
		min_uv_area = uv_area;
		max_uv_area = uv_area;
		no_cell = REF_FALSE;
	      }
	    else
	      {
		min_uv_area = MIN(min_uv_area,uv_area);
		max_uv_area = MAX(max_uv_area,uv_area);
	      }
	    if ( uv_area < 0.0 )
	      {
		n_neg++;
	      }
	    else
	      {
		n_pos++;
	      }
	  }
      if ( !no_cell )
	{
	  RSS( ref_geom_uv_area_sign( ref_grid, id, &sign_uv_area ), "sign");
	  printf ("face%5d: %4.1f %9.2e total (%10.3e,%10.3e) %d + %d -\n",
		  id, sign_uv_area, total_uv_area, min_uv_area, max_uv_area, n_pos, n_neg);
	}
    }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_inspect( REF_GEOM ref_geom )
{
  REF_INT geom;
  printf("ref_geom = %p\n",(void *)ref_geom);
  printf(" n = %d, max = %d\n",ref_geom_n(ref_geom),ref_geom_max(ref_geom));
  for ( geom = 0 ; geom < ref_geom_max( ref_geom ) ; geom++ )
    {
      switch(ref_geom_type(ref_geom,geom))
	{
	case REF_GEOM_NODE:
	  printf("%d node: %d global, %d id\n",
		 geom,
		 ref_geom_id(ref_geom,geom),
		 ref_geom_node(ref_geom,geom));
	  break;
	case REF_GEOM_EDGE:
	  printf("%d edge: %d id, %d global, t=%e\n",
		 geom,
		 ref_geom_id(ref_geom,geom),
		 ref_geom_node(ref_geom,geom),
		 ref_geom_param(ref_geom,0,geom));
	  break;
	case REF_GEOM_FACE:
	  printf("%d face: %d id, %d global, uv= %e %e\n",
		 geom,
		 ref_geom_id(ref_geom,geom),
		 ref_geom_node(ref_geom,geom),
		 ref_geom_param(ref_geom,0,geom),
		 ref_geom_param(ref_geom,1,geom) );
	  break;
	}
    }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tattle( REF_GEOM ref_geom, REF_INT node )
{
  REF_INT item, geom;

  printf(" tattle on node = %d\n",node);
  each_ref_adj_node_item_with_ref( ref_geom_adj(ref_geom), node, item, geom)
    {
      switch(ref_geom_type(ref_geom,geom))
	{
	case REF_GEOM_NODE:
	  printf("%d node: %d global, %d id\n",
		 geom,
		 ref_geom_id(ref_geom,geom),
		 ref_geom_node(ref_geom,geom));
	  break;
	case REF_GEOM_EDGE:
	  printf("%d edge: %d id, %d global, t=%e\n",
		 geom,
		 ref_geom_id(ref_geom,geom),
		 ref_geom_node(ref_geom,geom),
		 ref_geom_param(ref_geom,0,geom));
	  break;
	case REF_GEOM_FACE:
	  printf("%d face: %d id, %d global, uv= %e %e\n",
		 geom,
		 ref_geom_id(ref_geom,geom),
		 ref_geom_node(ref_geom,geom),
		 ref_geom_param(ref_geom,0,geom),
		 ref_geom_param(ref_geom,1,geom) );
	  break;
	}
    }
  
  return REF_SUCCESS;
}

REF_STATUS ref_geom_supported( REF_GEOM ref_geom, REF_INT node,
			       REF_BOOL *has_support )
{
  *has_support = !ref_adj_empty( ref_geom_adj(ref_geom), node );
  return REF_SUCCESS;
}

REF_STATUS ref_geom_add( REF_GEOM ref_geom, REF_INT node,
			 REF_INT type, REF_INT id,
			 REF_DBL *param )
{
  REF_INT geom;
  REF_INT orig, chunk;
  REF_INT max_limit = REF_INT_MAX/3;
  REF_STATUS status;

  if ( type < 0 || 2 < type )
    return REF_INVALID;

  status = ref_geom_find( ref_geom, node, type, id, &geom );
  RXS( status, REF_NOT_FOUND, "find failed");

  if ( REF_SUCCESS == status )
    {
      if ( type > 0 ) ref_geom_param(ref_geom,0,geom) = param[0];
      if ( type > 1 ) ref_geom_param(ref_geom,1,geom) = param[1];
      return REF_SUCCESS;
    }
    
  if ( REF_EMPTY == ref_geom_blank(ref_geom) )
    {
      RAS( ref_geom_max(ref_geom) != max_limit,
           "the number of geoms is too large for integers, cannot grow");
      orig = ref_geom_max(ref_geom);
      /* geometric growth for efficiency */
      chunk = MAX(1000,(REF_INT)( 1.5*(REF_DBL)orig ));

      /* try to keep under 32-bit limit */
      RAS( max_limit-orig > 0, "chunk limit at max");
      chunk = MIN( chunk, max_limit-orig );

      ref_geom_max(ref_geom) = orig + chunk;

      ref_realloc( ref_geom->descr, 3*ref_geom_max(ref_geom), REF_INT );
      ref_realloc( ref_geom->param, 2*ref_geom_max(ref_geom), REF_DBL );

      for (geom = orig; geom < ref_geom_max(ref_geom); geom++ )
        {
          ref_geom_type(ref_geom,geom) = REF_EMPTY;
          ref_geom_id(ref_geom,geom) = geom+1;
        }
      ref_geom_id(ref_geom,ref_geom_max(ref_geom)-1) = REF_EMPTY;
      ref_geom_blank(ref_geom) = orig;
    }

  geom = ref_geom_blank(ref_geom);
  ref_geom_blank(ref_geom) = ref_geom_id(ref_geom,geom);

  ref_geom_type(ref_geom,geom) = type;
  ref_geom_id(ref_geom,geom) = id;
  ref_geom_node(ref_geom,geom) = node;

  ref_geom_param(ref_geom,0,geom) = 0.0;
  ref_geom_param(ref_geom,1,geom) = 0.0;
  if ( type > 0 ) ref_geom_param(ref_geom,0,geom) = param[0];
  if ( type > 1 ) ref_geom_param(ref_geom,1,geom) = param[1];
  
  RSS( ref_adj_add(ref_geom->ref_adj, node, geom),"register geom" );

  ref_geom_n(ref_geom)++;

  return REF_SUCCESS;
}

REF_STATUS ref_geom_remove( REF_GEOM ref_geom, REF_INT node,
			    REF_INT type, REF_INT id)
{
  REF_INT geom;
  REF_STATUS status;
  
  status = ref_geom_find( ref_geom, node, type, id, &geom );
  RXS( status, REF_NOT_FOUND, "find failed");

  if ( REF_SUCCESS == status )
    {
      RSS( ref_adj_remove(ref_geom_adj(ref_geom),
			  ref_geom_node(ref_geom,geom), geom),
	   "unregister geom" );
	  
      ref_geom_type(ref_geom,geom) = REF_EMPTY;
      ref_geom_id(ref_geom,geom) = ref_geom_blank(ref_geom);
      ref_geom_blank(ref_geom) = geom;
      ref_geom_n(ref_geom)--;
    }
  
  return status;
}

REF_STATUS ref_geom_remove_all( REF_GEOM ref_geom, REF_INT node)
{
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_INT item, geom;

  item = ref_adj_first( ref_adj, node );
  while ( ref_adj_valid( item ) )
    {
      geom = ref_adj_item_ref( ref_adj, item );
      RSS( ref_adj_remove(ref_adj, node, geom),
	   "unregister geom" );
      
      ref_geom_type(ref_geom,geom) = REF_EMPTY;
      ref_geom_id(ref_geom,geom) = ref_geom_blank(ref_geom);
      ref_geom_blank(ref_geom) = geom;
      ref_geom_n(ref_geom)--;

      item = ref_adj_first( ref_adj, node );
    }
  
  return REF_SUCCESS;
}

REF_STATUS ref_geom_is_a( REF_GEOM ref_geom, REF_INT node,
			  REF_INT type, REF_BOOL *it_is )
{
  REF_INT item, geom;
  *it_is = REF_FALSE;
  each_ref_adj_node_item_with_ref( ref_geom_adj(ref_geom), node, item, geom)
    {
      if ( type == ref_geom_type(ref_geom,geom) )
	{
	  *it_is = REF_TRUE;
	  return REF_SUCCESS;
	}   
    }
  return REF_SUCCESS;
}

REF_STATUS ref_geom_unique_id( REF_GEOM ref_geom, REF_INT node,
			       REF_INT type, REF_INT *id)
{
  REF_INT item, geom;
  REF_BOOL found_one;
  found_one = REF_FALSE;
  each_ref_adj_node_item_with_ref( ref_geom_adj(ref_geom), node, item, geom)
    {
      if ( type == ref_geom_type(ref_geom,geom) )
	{
	  if ( found_one ) return REF_INVALID; /* second one makes invalid */
	  found_one = REF_TRUE;
	  *id = ref_geom_id(ref_geom,geom);
	}   
    }
  if ( found_one ) return REF_SUCCESS;
  return REF_NOT_FOUND;
}

REF_STATUS ref_geom_find( REF_GEOM ref_geom, REF_INT node,
			  REF_INT type, REF_INT id,
			  REF_INT *found )
{
  REF_INT item, geom;
  *found = REF_EMPTY;
  each_ref_adj_node_item_with_ref( ref_geom_adj(ref_geom), node, item, geom)
    {
      if ( type == ref_geom_type(ref_geom,geom) &&
	   id   == ref_geom_id(ref_geom,geom) )
	{
	  *found = geom;
	  return REF_SUCCESS;
	}   
    }
  return REF_NOT_FOUND;
}

REF_STATUS ref_geom_tuv( REF_GEOM ref_geom, REF_INT node,
			 REF_INT type, REF_INT id,
			 REF_DBL *param )
{
  REF_INT geom;
  
  RSS( ref_geom_find( ref_geom, node, type, id, &geom ), "not found");

  if ( type > 0 ) param[0] = ref_geom_param(ref_geom,0,geom);
  if ( type > 1 ) param[1] = ref_geom_param(ref_geom,1,geom);

  return REF_SUCCESS;
}

REF_STATUS ref_geom_add_between( REF_GRID ref_grid,
				 REF_INT node0, REF_INT node1, 
				 REF_INT new_node )
{
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT item0, item1;
  REF_INT geom0, geom1;
  REF_INT type, id;
  REF_DBL param[2], param0[2], param1[2];
  REF_BOOL has_id;

  each_ref_adj_node_item_with_ref( ref_geom_adj(ref_geom),
				   node0, item0, geom0)
    each_ref_adj_node_item_with_ref( ref_geom_adj(ref_geom),
				     node1, item1, geom1)
    {
      type = REF_GEOM_EDGE;
      if ( ref_geom_type(ref_geom,geom0) == type &&
	   ref_geom_type(ref_geom,geom1) == type &&
	   ref_geom_id(ref_geom,geom0) == ref_geom_id(ref_geom,geom1) )
	{
	  id = ref_geom_id(ref_geom,geom0);
	  RSS(ref_cell_side_has_id(ref_grid_edg(ref_grid),node0,node1,id,
				   &has_id), "has edge id");
	  if ( has_id )
	    {
	      RSS( ref_geom_tuv(ref_geom,node0,type,id,param0), "node0" );
	      RSS( ref_geom_tuv(ref_geom,node1,type,id,param1), "node1" );
	      if ( type > 0 ) param[0] = 0.5 * ( param0[0] + param1[0] );
	      if ( type > 1 ) param[1] = 0.5 * ( param0[1] + param1[1] );
	      RSS( ref_geom_add(ref_geom,new_node,type,id,param), "new geom" );
	    }
	}
      type = REF_GEOM_FACE;
      if ( ref_geom_type(ref_geom,geom0) == type &&
	   ref_geom_type(ref_geom,geom1) == type &&
	   ref_geom_id(ref_geom,geom0) == ref_geom_id(ref_geom,geom1) )
	{
	  id = ref_geom_id(ref_geom,geom0);
	  RSS(ref_cell_side_has_id(ref_grid_tri(ref_grid),node0,node1,id,
				   &has_id), "has edge id");
	  if ( has_id )
	    {
	      RSS( ref_geom_tuv(ref_geom,node0,type,id,param0), "node0" );
	      RSS( ref_geom_tuv(ref_geom,node1,type,id,param1), "node1" );
	      if ( type > 0 ) param[0] = 0.5 * ( param0[0] + param1[0] );
	      if ( type > 1 ) param[1] = 0.5 * ( param0[1] + param1[1] );
	      RSS( ref_geom_add(ref_geom,new_node,type,id,param), "new geom" );
	    }
	}
    }
    
  return REF_SUCCESS;
}

REF_STATUS ref_geom_eval_edge_face_uv( REF_GEOM ref_geom, REF_INT edge_geom )
{
#ifdef HAVE_EGADS
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_INT node, item, face_geom;
  double t;
  double uv[2];
  int sense = 0;
  ego *edges, *faces;
  ego edge, face;

  if ( edge_geom < 0 || ref_geom_max(ref_geom) <= edge_geom )
    return REF_INVALID;
  if ( REF_GEOM_EDGE != ref_geom_type(ref_geom,edge_geom) )
    return REF_INVALID;

  edges = (ego *)(ref_geom->edges);
  edge = edges[ref_geom_id(ref_geom,edge_geom) - 1]; 

  t = ref_geom_param(ref_geom,0,edge_geom);

  node = ref_geom_node(ref_geom,edge_geom);

  faces = (ego *)(ref_geom->faces);
  each_ref_adj_node_item_with_ref( ref_adj, node, item, face_geom)
    {
      if (REF_GEOM_FACE == ref_geom_type(ref_geom,face_geom))
	{
	  face = faces[ref_geom_id(ref_geom,face_geom) - 1];
	  REIS( EGADS_SUCCESS,
		EG_getEdgeUV(face, edge, sense, t, uv), "eval edge face uv");
	  ref_geom_param(ref_geom,0,face_geom) = uv[0];
	  ref_geom_param(ref_geom,1,face_geom) = uv[1];
	}
    }

  return REF_SUCCESS;
#else
  if ( edge_geom < 0 || ref_geom_max(ref_geom) <= edge_geom )
    return REF_INVALID;
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_geom_constrain( REF_GRID ref_grid, REF_INT node )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_INT item, geom;
  REF_BOOL have_geom_node;
  REF_BOOL have_geom_edge;
  REF_BOOL have_geom_face;
  REF_INT edge_geom;
  REF_INT face_geom;
  REF_DBL xyz[3];

  /* put sloppy geom handling here */
  
  /* no geom, do nothing */
  if (ref_adj_empty( ref_adj, node )) return REF_SUCCESS;

  have_geom_node = REF_FALSE;
  each_ref_adj_node_item_with_ref( ref_adj, node, item, geom)
    {
      if (REF_GEOM_NODE == ref_geom_type(ref_geom,geom))
	{
	  have_geom_node = REF_TRUE;
	  break;
	}
    }

  /* node geom, do nothing */
  if (have_geom_node)
    {
      return REF_SUCCESS;
    }
  
  have_geom_edge = REF_FALSE;
  edge_geom = REF_EMPTY;
  each_ref_adj_node_item_with_ref( ref_adj, node, item, geom)
    {
      if (REF_GEOM_EDGE == ref_geom_type(ref_geom,geom))
	{
	  have_geom_edge = REF_TRUE;
	  edge_geom = geom;
	  break;
	}
    }
  
  /* edge geom, evaluate edge and update face uv */
  if (have_geom_edge)
    {
      RSS( ref_geom_eval( ref_geom, edge_geom, xyz, NULL ), "eval edge" );
      node = ref_geom_node(ref_geom,edge_geom);
      ref_node_xyz(ref_node,0,node) = xyz[0];
      ref_node_xyz(ref_node,1,node) = xyz[1];
      ref_node_xyz(ref_node,2,node) = xyz[2];
      RSS( ref_geom_eval_edge_face_uv( ref_geom, edge_geom ), "resol edge uv");
      return REF_SUCCESS;
    }
  
  /* look for face geom */
  have_geom_face = REF_FALSE;
  face_geom = REF_EMPTY;
  each_ref_adj_node_item_with_ref( ref_adj, node, item, geom)
    {
      if (REF_GEOM_FACE == ref_geom_type(ref_geom,geom))
	{
	  have_geom_face = REF_TRUE;
	  face_geom = geom;
	  break;
	}
    }

  /* face geom, evaluate on face uv */
  if (have_geom_face)
    {
      RSS( ref_geom_eval( ref_geom, face_geom, xyz, NULL ), "eval face" );
      node = ref_geom_node(ref_geom,face_geom);
      ref_node_xyz(ref_node,0,node) = xyz[0];
      ref_node_xyz(ref_node,1,node) = xyz[1];
      ref_node_xyz(ref_node,2,node) = xyz[2];
      return REF_SUCCESS;
    }

  return REF_SUCCESS;
}
/* 
   [x_t,y_t,z_t] edge
   [x_tt,y_tt,z_tt]
   [x_u,y_u,z_u] [x_v,y_v,z_v] face
   [x_uu,y_uu,z_uu] [x_uv,y_uv,z_uv] [x_vv,y_vv,z_vv]
*/
REF_STATUS ref_geom_eval( REF_GEOM ref_geom, REF_INT geom,
			  REF_DBL *xyz, REF_DBL *dxyz_dtuv )
{
#ifdef HAVE_EGADS
  double eval[18];
  double params[2];
  REF_INT i;
  ego *nodes, *edges, *faces;
  ego object;
  if ( geom < 0 || ref_geom_max(ref_geom) <= geom )
    return REF_INVALID;
  params[0] = 0.0; params[1] = 0.0;
  object = (ego)NULL;
  switch (ref_geom_type(ref_geom,geom))
    {
    case (REF_GEOM_NODE) :
      RNS(ref_geom->nodes,"nodes not loaded");
      if ( ref_geom_id(ref_geom,geom) < 1 ||
	   ref_geom_id(ref_geom,geom) > ref_geom->nnode ) 
	return REF_INVALID;
      nodes = (ego *)(ref_geom->nodes);
      object = nodes[ref_geom_id(ref_geom,geom) - 1]; 
      {
	ego ref, *pchldrn;
	int oclass, mtype, nchild, *psens;
	REIS( EGADS_SUCCESS,
	      EG_getTopology(object, &ref, &oclass, &mtype, xyz,
			     &nchild, &pchldrn, &psens), "EG topo node");
      }
      return REF_SUCCESS;
      break;
    case (REF_GEOM_EDGE) :
      RNS(ref_geom->edges,"edges not loaded");
      if ( ref_geom_id(ref_geom,geom) < 1 ||
	   ref_geom_id(ref_geom,geom) > ref_geom->nedge ) 
	return REF_INVALID;
      edges = (ego *)(ref_geom->edges);
      object = edges[ref_geom_id(ref_geom,geom) - 1]; 
      params[0] = ref_geom_param(ref_geom,0,geom);
      break;
    case (REF_GEOM_FACE) :
      RNS(ref_geom->faces,"faces not loaded");
      if ( ref_geom_id(ref_geom,geom) < 1 ||
	   ref_geom_id(ref_geom,geom) > ref_geom->nface ) 
	return REF_INVALID;
      faces = (ego *)(ref_geom->faces);
      object = faces[ref_geom_id(ref_geom,geom) - 1]; 
      params[0] = ref_geom_param(ref_geom,0,geom);
      params[1] = ref_geom_param(ref_geom,1,geom);
      break;
    default:
      RSS(REF_IMPLEMENT, "unknown geom" );
    }
  
  REIS( EGADS_SUCCESS,
	EG_evaluate(object, params, eval), "eval");
  xyz[0]=eval[0];
  xyz[1]=eval[1];
  xyz[2]=eval[2];
  if ( NULL != dxyz_dtuv )
    {
      for (i=0;i<6;i++) dxyz_dtuv[i] = eval[3+i];
      if (REF_GEOM_FACE == ref_geom_type(ref_geom,geom))
	for (i=0;i<9;i++) dxyz_dtuv[6+i] = eval[9+i];
    }
  return REF_SUCCESS;
#else
  if ( geom < 0 || ref_geom_max(ref_geom) <= geom )
    return REF_INVALID;
  printf("evaluating to (0,0,0), No EGADS linked for %s\n", __func__);
  xyz[0] = 0.0;
  xyz[1] = 0.0;
  xyz[2] = 0.0;
  if ( NULL != dxyz_dtuv ) dxyz_dtuv[0] = 0.0;
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_geom_curvature( REF_GEOM ref_geom, REF_INT geom,
			       REF_DBL *kr, REF_DBL *r,
			       REF_DBL *ks, REF_DBL *s )
{
#ifdef HAVE_EGADS
  double curvature[8];
  double params[2];
  ego *faces;
  ego object;
  int egads_status;
  if ( geom < 0 || ref_geom_max(ref_geom) <= geom )
    return REF_INVALID;
  params[0] = 0.0; params[1] = 0.0;
  object = (ego)NULL;
  switch (ref_geom_type(ref_geom,geom))
    {
    case (REF_GEOM_NODE) :
      RSS(REF_IMPLEMENT, "geom node" );
      break;
    case (REF_GEOM_EDGE) :
      RSS(REF_IMPLEMENT, "geom edge" );
      break;
    case (REF_GEOM_FACE) :
      RNS(ref_geom->faces,"faces not loaded");
      faces = (ego *)(ref_geom->faces);
      object = faces[ref_geom_id(ref_geom,geom) - 1]; 
      params[0] = ref_geom_param(ref_geom,0,geom);
      params[1] = ref_geom_param(ref_geom,1,geom);
      break;
    default:
      RSS(REF_IMPLEMENT, "unknown geom" );
    }
  RNS(object,"EGADS object is NULL. Has the geometry been loaded?");
  
  egads_status = EG_curvature(object, params, curvature);
  if ( EGADS_DEGEN == egads_status )
    {
      REF_DBL xyz[3], dxyz_duv[15], du, dv;
      ego ref, *pchldrn;
      int oclass, mtype, nchild, *psens;
      double uv_range[4];
      REF_DBL shift=1.0e-4;
      RSS( ref_geom_eval( ref_geom, geom, xyz, dxyz_duv ), "eval" );
      du = sqrt(ref_math_dot(&(dxyz_duv[0]),&(dxyz_duv[0])));
      dv = sqrt(ref_math_dot(&(dxyz_duv[3]),&(dxyz_duv[3])));
      REIS( EGADS_SUCCESS,
	    EG_getTopology(object, &ref, &oclass, &mtype, uv_range,
			   &nchild, &pchldrn, &psens), "EG topo face");
      if (du > dv)
	{
	  params[0] = (1.0-shift)*params[0]+shift*0.5*(uv_range[0]+uv_range[1]);
	}
      else
	{
	  params[1] = (1.0-shift)*params[1]+shift*0.5*(uv_range[2]+uv_range[3]);
	}
      egads_status = EG_curvature(object, params, curvature);
    }
  REIS( EGADS_SUCCESS, egads_status, "curve");
  *kr=curvature[0];
  r[0] = curvature[1];
  r[1] = curvature[2];
  r[2] = curvature[3];
  *ks=curvature[4];
  s[0] = curvature[5];
  s[1] = curvature[6];
  s[2] = curvature[7];
  return REF_SUCCESS;
#else
  if ( geom < 0 || ref_geom_max(ref_geom) <= geom )
    return REF_INVALID;
  printf("curvature 0, 0: No EGADS linked for %s\n", __func__);
  *kr=0.0;
  r[0] = 1.0;
  r[1] = 0.0;
  r[2] = 0.0;
  *ks=0.0;
  s[0] = 0.0;
  s[1] = 1.0;
  s[2] = 0.0;
  return REF_IMPLEMENT;
#endif
}

REF_STATUS ref_geom_uv_rsn( REF_DBL *uv,
			    REF_DBL *r, REF_DBL *s, REF_DBL *n,
			    REF_DBL *drsduv )
{
  REF_INT i;
  REF_DBL dot;
  REF_DBL len;

  for (i=0;i<3;i++) r[i] = uv[i];
  drsduv[0]=1.0;drsduv[1]=0.0;
  for (i=0;i<3;i++) s[i] = uv[i+3];
  drsduv[2]=0.0;drsduv[3]=1.0;
  len = sqrt(ref_math_dot(r,r));
  drsduv[0] /= len;
  drsduv[1] /= len;
  RSS( ref_math_normalize( r ), "norm r (u)" );
  len = sqrt(ref_math_dot(s,s));
  drsduv[2] /= len;
  drsduv[3] /= len;  
  RSS( ref_math_normalize( s ), "norm s (v)" );

  dot = ref_math_dot(r,s);
  for (i=0;i<3;i++) s[i] -= dot*r[i];
  drsduv[2] -= dot*drsduv[0];
  drsduv[3] -= dot*drsduv[1];
  
  len = sqrt(ref_math_dot(s,s));
  drsduv[2] /= len;
  drsduv[3] /= len;  
  RSS( ref_math_normalize( s ), "norm s (v)" );

  ref_math_cross_product( r, s, n );

  return REF_SUCCESS;
}

REF_STATUS ref_geom_rsn( REF_GEOM ref_geom, REF_INT geom,
			 REF_DBL *r, REF_DBL *s, REF_DBL *n )
{
  REF_DBL kr,ks;
  RSS( ref_geom_curvature( ref_geom, geom, &kr, r, &ks, s ), "eval face" );
  ref_math_cross_product( r, s, n );
  return REF_SUCCESS;
}

REF_STATUS ref_geom_verify_param( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT geom;
  REF_INT node;
  REF_DBL xyz[3];
  REF_DBL dist;
  REF_DBL tol;
  REF_BOOL geom_edge;
  
  each_ref_geom_node( ref_geom, geom )
    {
      RSS( ref_geom_eval( ref_geom, geom, xyz, NULL ), "eval xyz" );
      node = ref_geom_node(ref_geom,geom);
      dist = sqrt( pow(xyz[0]-ref_node_xyz(ref_node,0,node),2) +
		   pow(xyz[1]-ref_node_xyz(ref_node,1,node),2) +
		   pow(xyz[2]-ref_node_xyz(ref_node,2,node),2) );
      tol = 1.0e-12;
      if ( dist > tol )
	{
	  printf("geom node %d node %d dist %e\n",geom,node,dist);	 
	  RSS( ref_geom_tattle( ref_geom, node ), "tattle");
	}
    }
  
  each_ref_geom_edge( ref_geom, geom )
    {
      RSS( ref_geom_eval( ref_geom, geom, xyz, NULL ), "eval xyz" );
      node = ref_geom_node(ref_geom,geom);
      dist = sqrt( pow(xyz[0]-ref_node_xyz(ref_node,0,node),2) +
		   pow(xyz[1]-ref_node_xyz(ref_node,1,node),2) +
		   pow(xyz[2]-ref_node_xyz(ref_node,2,node),2) );
      tol = 1.0e-12;
      if ( dist > tol )
	{
	  printf("geom edge %d node %d dist %e\n",geom,node,dist);	 
	  RSS( ref_geom_tattle( ref_geom, node ), "tattle");
	}
    }
  
  each_ref_geom_face( ref_geom, geom )
    {
      RSS( ref_geom_eval( ref_geom, geom, xyz, NULL ), "eval xyz" );
      node = ref_geom_node(ref_geom,geom);
      RSS( ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &geom_edge), "edge");
      dist = sqrt( pow(xyz[0]-ref_node_xyz(ref_node,0,node),2) +
		   pow(xyz[1]-ref_node_xyz(ref_node,1,node),2) +
		   pow(xyz[2]-ref_node_xyz(ref_node,2,node),2) );
      tol = 1.0e-12;
      if (geom_edge) tol = 2.0e-5;
      if ( dist > tol )
	{
	  printf("geom face %d node %d dist %e\n",geom,node,dist);	 
	  RSS( ref_geom_tattle( ref_geom, node ), "tattle");
	}
    }
  
  return REF_SUCCESS;
}

REF_STATUS ref_geom_verify_topo( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT node;
  REF_BOOL geom_node, geom_edge, geom_face;
  REF_BOOL no_face, no_edge;
  
  for ( node = 0 ;node < ref_node_max(ref_node);node++ )
    if ( ref_node_valid( ref_node, node ) )
      {
	RSS( ref_geom_is_a(ref_geom, node, REF_GEOM_NODE, &geom_node), "node");
	RSS( ref_geom_is_a(ref_geom, node, REF_GEOM_EDGE, &geom_edge), "edge");
	RSS( ref_geom_is_a(ref_geom, node, REF_GEOM_FACE, &geom_face), "face");
	no_face =
	  ref_cell_node_empty( ref_grid_tri( ref_grid ), node ) &&
	  ref_cell_node_empty( ref_grid_qua( ref_grid ), node );
	no_edge =
	  ref_cell_node_empty( ref_grid_edg( ref_grid ), node );
	if ( geom_node )
	  {
	    if (no_edge && ref_mpi_id == ref_node_part(ref_node,node) ) 
	      {
		THROW("geom node missing edge");
	      }
	    if (no_face && ref_mpi_id == ref_node_part(ref_node,node) ) {
	      THROW("geom node missing tri or qua");
	    }
	}
	if ( geom_edge )
	  {
	    if (no_edge && ref_mpi_id == ref_node_part(ref_node,node) ) 
	      {
		RSS(ref_node_location(ref_node,node),"loc");
		RSS(ref_geom_tattle(ref_geom,node),"tatt");
		RSS(ref_geom_tec( ref_grid, "ref_geom_typo_error.tec" ),
		    "geom tec" );
		THROW("geom edge missing edge");
	      }
	    if (no_face && ref_mpi_id == ref_node_part(ref_node,node) )
	      { 
		RSS(ref_node_location(ref_node,node),"loc");
		RSS(ref_geom_tattle(ref_geom,node),"tatt");
		RSS(ref_geom_tec( ref_grid, "ref_geom_typo_error.tec" ),
		    "geom tec" );
		THROW("geom edge missing tri or qua");
	      }
	  }
	if ( geom_face )
	  {
	    if (no_face && ref_mpi_id == ref_node_part(ref_node,node) )
	      {
		printf("no face for geom\n");
		RSS(ref_node_location(ref_node,node),"loc");
		RSS(ref_geom_tattle(ref_geom,node),"tatt");
		RSS(ref_geom_tec( ref_grid, "ref_geom_typo_error.tec" ),
		    "geom tec" );
		THROW("geom face missing tri or qua");
	      }
	  }
	if ( !no_face )
	  {
	    if ( !geom_face)
	      {
		printf("no geom for face\n");
		RSS(ref_node_location(ref_node,node),"loc");
		RSS(ref_geom_tattle(ref_geom,node),"tatt");
		RSS(ref_geom_tec( ref_grid, "ref_geom_typo_error.tec" ),
		    "geom tec" );
		THROW("geom face missing tri or qua");
	      }
	  }
	if ( geom_face && !geom_edge )
	  {
	    REF_INT item, geom;
	    REF_BOOL found_one;
	    REF_BOOL found_too_many;
	    found_one = REF_FALSE;
	    found_too_many = REF_FALSE;
	    each_ref_adj_node_item_with_ref( ref_geom_adj(ref_geom),
					     node, item, geom)
	      {
		if ( REF_GEOM_FACE == ref_geom_type(ref_geom,geom) )
		  {
		    if ( found_one ) found_too_many = REF_TRUE;
		    found_one = REF_TRUE;
		  }   
	      }
	    if (!found_one || found_too_many)
	      {
		if (!found_one) printf("none found\n");
		if (found_too_many) printf("found too many\n");
		RSS(ref_node_location(ref_node,node),"loc");
		RSS(ref_geom_tattle(ref_geom,node),"tatt");
		RSS(ref_geom_tec( ref_grid, "ref_geom_typo_error.tec" ),
		    "geom tec" );
		THROW("geom face missing tri or qua");
	      }
	  }
      }
    else
      {
	if (!ref_adj_empty( ref_geom_adj(ref_geom), node ) )
	  THROW("invalid node has geom");
      }
     
    return REF_SUCCESS;
}

REF_STATUS ref_geom_tetgen_volume( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  char *smesh_name = "ref_geom_test.smesh";
  char *node_name = "ref_geom_test.1.node";
  char *ele_name = "ref_geom_test.1.ele";
  char command[1024];
  FILE *file;
  REF_INT nnode, ndim, attr, mark;
  REF_INT ntet, node_per;
  REF_INT node, nnode_surface, item, new_node;
  REF_DBL xyz[3], dist;
  REF_INT cell, new_cell, nodes[REF_CELL_MAX_SIZE_PER];
  int system_status;
  RSS( ref_export_smesh( ref_grid, smesh_name ), "smesh" );
  sprintf( command, "tetgen -pYq1.0/0z %s > %s.out", smesh_name, smesh_name );
  system_status = system( command );
  if ( 0 != system_status )
    {
      printf("tec360 ref_geom_test_debug_geom.tec\n");
      RSS(ref_geom_tec( ref_grid, "ref_geom_test_debug_geom.tec" ),
	  "dbg geom" );
      printf("tec360 ref_geom_test_debug_surf.tec\n");
      RSS( ref_export_tec_surf( ref_grid,
				"ref_geom_test_debug_surf.tec"), "dbg surf" );
    }
  REIS(0, system_status, "tetgen failed");

  file = fopen(node_name,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",node_name);
  RNS(file, "unable to open file" );

  REIS( 1, fscanf( file, "%d", &nnode ), "node header nnode" );
  REIS( 1, fscanf( file, "%d", &ndim ), "node header ndim" );
  REIS( 3, ndim, "not 3D");
  REIS( 1, fscanf( file, "%d", &attr ), "node header attr" );
  REIS( 0, attr, "nodes have attribute 3D");
  REIS( 1, fscanf( file, "%d", &mark ), "node header mark" );
  REIS( 0, mark, "nodes have mark");

  /* verify surface nodes */
  nnode_surface = ref_node_n(ref_node);
  for( node=0; node<nnode_surface ; node++ ) 
    {
      REIS( 1, fscanf( file, "%d", &item ), "node item" );
      RES( node, item, "node index");
      RES( 1, fscanf( file, "%lf", &(xyz[0]) ), "x" );
      RES( 1, fscanf( file, "%lf", &(xyz[1]) ), "y" );
      RES( 1, fscanf( file, "%lf", &(xyz[2]) ), "z" );
      dist = sqrt( (xyz[0]-ref_node_xyz( ref_node, 0, node )) *
		   (xyz[0]-ref_node_xyz( ref_node, 0, node )) +
		   (xyz[1]-ref_node_xyz( ref_node, 1, node )) *
		   (xyz[1]-ref_node_xyz( ref_node, 1, node )) +
		   (xyz[2]-ref_node_xyz( ref_node, 2, node )) *
		   (xyz[2]-ref_node_xyz( ref_node, 2, node )) );
      if ( dist > 1.0e-12 )
	{
	  printf("node %d off by %e\n",node,dist);
	  THROW("tetgen moved node");
	}
    }

  /* interior nodes */
  for( node=nnode_surface; node<nnode ; node++ ) 
    {
      REIS( 1, fscanf( file, "%d", &item ), "node item" );
      RES( node, item, "node index");
      RSS( ref_node_add(ref_node, node, &new_node ), "new_node");
      RES( node, new_node, "node index");
      RES( 1, fscanf( file, "%lf", &(xyz[0]) ), "x" );
      RES( 1, fscanf( file, "%lf", &(xyz[1]) ), "y" );
      RES( 1, fscanf( file, "%lf", &(xyz[2]) ), "z" );
      ref_node_xyz( ref_node, 0, new_node ) = xyz[0];
      ref_node_xyz( ref_node, 1, new_node ) = xyz[1];
      ref_node_xyz( ref_node, 2, new_node ) = xyz[2];
    }
  
  fclose( file );
  
  /* check faces when paranoid, but tetgen -z should not mess with them */

  file = fopen(ele_name,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",ele_name);
  RNS(file, "unable to open file" );

  REIS( 1, fscanf( file, "%d", &ntet ), "ele header ntet" );
  REIS( 1, fscanf( file, "%d", &node_per ), "ele header node_per" );
  REIS( 4, node_per, "expected tets");
  REIS( 1, fscanf( file, "%d", &mark ), "ele header mark" );
  REIS( 0, mark, "ele have mark");

  ref_cell = ref_grid_tet(ref_grid);
  for( cell = 0; cell < ntet ; cell++ )
    {
      REIS( 1, fscanf( file, "%d", &item ), "tet item" );
      RES( cell, item, "node index");
      for ( node = 0 ; node < 4 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "tet" );
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "new tet");
      RES( cell, new_cell, "tet index");
    }
  
  fclose( file );
  
  return REF_SUCCESS;
}
  
REF_STATUS ref_geom_egads_load( REF_GEOM ref_geom, const char *filename )
{
#ifdef HAVE_EGADS
  ego context;
  ego model = NULL;
  ego geom, *bodies, *children;
  int oclass, mtype, nbody, *senses, nchild;
  ego solid, *faces, *edges, *nodes;
  int nface, nedge, nnode;

  context = (ego)(ref_geom->context);

  REIS( EGADS_SUCCESS, EG_loadModel(context, 0, filename, &model), "EG load");

  /*

#ifndef HAVE_EGADS_LITE

  REIS( EGADS_SUCCESS, EG_loadModel(context, 0, filename, &model), "EG load");

  {
    long cad_data_size;
    REF_BYTE *cad_data;
    int EG_exportModel(ego model, long *cad_data_size, char **cad_data);
    REIS( EGADS_SUCCESS, EG_exportModel(model, &cad_data_size, &cad_data), 
	  "EG stream");
    ref_geom_cad_data_size(ref_geom) = cad_data_size;
    ref_malloc(ref_geom_cad_data(ref_geom), ref_geom_cad_data_size(ref_geom),
	       REF_BYTE );
    memcpy( ref_geom_cad_data(ref_geom), cad_data, 
	    ref_geom_cad_data_size(ref_geom) );
    free( cad_data );
  }


#else
  if (ref_mpi_master )
    printf("EGADS lite, using meshb data ignore %s\n",filename);

  RAS( 0 < ref_geom_cad_data_size(ref_geom), "zero size cad_data" );
  RNS( ref_geom_cad_data(ref_geom), "cad_data NULL" );
  REIS( EGADS_SUCCESS, EG_importModel(context, 
				      (long)ref_geom_cad_data_size(ref_geom),
				      ref_geom_cad_data(ref_geom),
				      &model ), "EG load");
#endif

  */
 
  REIS( EGADS_SUCCESS,
	EG_getTopology(model, &geom, &oclass, &mtype, NULL,
		       &nbody, &bodies, &senses), "EG topo bodies");
  REIS( 1, nbody, "expected 1 body" );
  solid = bodies[0];
  REIS( EGADS_SUCCESS,
	EG_getTopology(solid, &geom, &oclass, &mtype,
		       NULL, &nchild, &children, &senses), "EG topo body type");
  REIS( SOLIDBODY, mtype, "expected SOLIDBODY" );
  ref_geom->solid = (void *)solid;
  
  REIS( EGADS_SUCCESS,
	EG_getBodyTopos(solid, NULL, NODE, &nnode, &nodes), "EG node topo");
  ref_geom->nnode = nnode;
  ref_geom->nodes = (void *)nodes;
  REIS( EGADS_SUCCESS,
	EG_getBodyTopos(solid, NULL, EDGE, &nedge, &edges), "EG edge topo");
  ref_geom->nedge = nedge;
  ref_geom->edges = (void *)edges;
  REIS( EGADS_SUCCESS,
	EG_getBodyTopos(solid, NULL, FACE, &nface, &faces), "EG face topo");
  ref_geom->nface = nface;
  ref_geom->faces = (void *)faces;
  
#else
  printf("returning empty grid from %s, No EGADS linked for %s\n",
	 __func__,filename);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
#endif
  
  return REF_SUCCESS;
}

REF_STATUS ref_geom_infer_nedge_nface( REF_GRID ref_grid )
{
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT min_id, max_id;
  RSS( ref_geom_faceid_range( ref_grid, &min_id, &max_id), "face range");
  REIS( 1, min_id, "first face id not 1" );
  ref_geom->nface = max_id;
  RSS( ref_geom_edgeid_range( ref_grid, &min_id, &max_id), "edge range");
  REIS( 1, min_id, "first edge id not 1" );
  ref_geom->nedge = max_id;
  return REF_SUCCESS;
}

REF_STATUS ref_geom_egads_diagonal( REF_GEOM ref_geom, REF_DBL *diag )
{
#ifdef HAVE_EGADS
  ego solid;
  double box[6];
  solid = (ego)(ref_geom->solid);
  
  RNS(solid,"EGADS solid object is NULL. Has the geometry been loaded?");
  
  REIS( EGADS_SUCCESS, EG_getBoundingBox(solid, box), "EG bounding box");
  *diag = sqrt((box[0]-box[3])*(box[0]-box[3]) +
	       (box[1]-box[4])*(box[1]-box[4]) +
	       (box[2]-box[5])*(box[2]-box[5]));

#else
  printf("returning 1.0 from %s, No EGADS\n",
	 __func__);
  *diag = 1.0;
  SUPRESS_UNUSED_COMPILER_WARNING(ref_geom);
#endif
  
  return REF_SUCCESS;
}

REF_STATUS ref_geom_egads_tess( REF_GRID ref_grid, REF_DBL max_length )
{
#ifdef HAVE_EGADS
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tri, new_cell;
  REF_DBL param[2];
  double params[3], size;
  ego geom;
  ego solid, tess;
  int tess_status, nvert;
  int face, edge, plen, tlen;
  const double *points, *uv, *t;
  const int    *ptype, *pindex, *tris, *tric;
  int node, new_node, pty, pin;
  double verts[3];
  
  solid = (ego)(ref_geom->solid);

  RSS( ref_geom_egads_diagonal( ref_geom, &size ), "bbox diag");
  /* maximum length of an EDGE segment or triangle side (in physical space) */
  if ( max_length > 0.0 )
    {
      params[0] = max_length;
    }
  else
    {
      params[0] = ABS(max_length)*size;
    }
  /* curvature-based value that looks locally at the deviation between
     the centroid of the discrete object and the underlying geometry */
  params[1] =  0.001*size;
  /* maximum interior dihedral angle (in degrees) */
  params[2] = 15.0;
  /* printf("params = %f,%f,%f\n",params[0],params[1],params[2]); */

  REIS( EGADS_SUCCESS,
	EG_makeTessBody(solid, params, &tess), "EG tess");
  REIS( EGADS_SUCCESS,
	EG_statusTessBody(tess, &geom, &tess_status, &nvert), "EG tess");
  REIS( 1, tess_status, "tess not closed" );

  for (node = 0; node < nvert; node++) {
    REIS( EGADS_SUCCESS,
	  EG_getGlobal(tess, node+1, &pty, &pin, verts), "global node info");
    RSS( ref_node_add(ref_node, node, &new_node ), "new_node");
    RES( node, new_node, "node index");
    ref_node_xyz( ref_node, 0, node ) = verts[0];
    ref_node_xyz( ref_node, 1, node ) = verts[1];
    ref_node_xyz( ref_node, 2, node ) = verts[2];
    if ( 0 == pty )
      RSS( ref_geom_add( ref_geom, node, REF_GEOM_NODE, pin, NULL), "node");
  }

  for (face = 0; face < (ref_geom->nface); face++) {
    REIS( EGADS_SUCCESS,
	  EG_getTessFace(tess, face+1, &plen, &points, &uv, &ptype, &pindex,
			 &tlen, &tris, &tric), "tess query face" );
    
    for ( node = 0; node<plen; node++ ) {
      REIS( EGADS_SUCCESS,
	    EG_localToGlobal(tess, face+1, node+1, &(nodes[0])), "l2g0");
      nodes[0] -= 1;
      param[0] = uv[0+2*node];
      param[1] = uv[1+2*node];
      RSS( ref_geom_add( ref_geom, nodes[0], REF_GEOM_FACE, face+1, param),
	   "face uv");
    }
    for ( tri = 0; tri<tlen; tri++ ) {
      REIS( EGADS_SUCCESS,
	    EG_localToGlobal(tess, face+1, tris[0+3*tri], &(nodes[0])), "l2g0");
      REIS( EGADS_SUCCESS,
	    EG_localToGlobal(tess, face+1, tris[1+3*tri], &(nodes[1])), "l2g1");
      REIS( EGADS_SUCCESS,
	    EG_localToGlobal(tess, face+1, tris[2+3*tri], &(nodes[2])), "l2g2");
      nodes[0] -= 1;
      nodes[1] -= 1;
      nodes[2] -= 1;
      nodes[3] = face + 1;
      RSS( ref_cell_add(ref_grid_tri(ref_grid), nodes, &new_cell ), "new tri");
    }


    {
      ego esurf, *eloops, eref;
      int oclass, mtype, nloop,*senses,*pinfo;
      double trange[2], data[18], *preal;
      REIS( EGADS_SUCCESS,
	    EG_getTopology(((ego *)(ref_geom->faces))[face],
			   &esurf, &oclass, &mtype,
			   trange, &nloop, &eloops, &senses), "topo" );
      REIS( EGADS_SUCCESS,
	    EG_getGeometry(esurf, &oclass, &mtype,
			   &eref, &pinfo, &preal),"geom");
      EG_free(pinfo);
      EG_free(preal);
      if (mtype != PLANE)
	{
	  ego ecurve, *eedges, *echilds;
	  int iloop, iedge, nedge, nchild, inode;
	  /* loop through all Loops associated with this Face */
	  for (iloop = 0; iloop < nloop; iloop++)
	    {
	      /* loop through all Edges associated with this Loop */
	      REIS( EGADS_SUCCESS,
		    EG_getTopology(eloops[iloop], &ecurve, &oclass, &mtype,
				   data, &nedge, &eedges, &senses), "topo");
	      for (iedge = 0; iedge < nedge; iedge++)
		{
		  REIS( EGADS_SUCCESS,
			EG_getTopology(eedges[iedge], &ecurve, &oclass, &mtype,
				       trange, &nchild, &echilds, &senses), "tp");
		  if (mtype == DEGENERATE)
		    {
		      double uvmin[6], uvmax[6];
		      printf("face id %d has degen\n",face+1);
		      /* find index of bounding Node */
		      inode = EG_indexBodyTopo(solid, echilds[0]);
		      REIS( EGADS_SUCCESS,
			    EG_evaluate(eedges[iedge+nedge], &trange[0], uvmin),
			    "eval min");
		      REIS( EGADS_SUCCESS,
			    EG_evaluate(eedges[iedge+nedge], &trange[1], uvmax),
			    "eval max");
		      for (node = 0; node<plen; node++)
			{
			  if (ptype[node] == 0 && pindex[node] == inode)
			    {
			      REIS( EGADS_SUCCESS,
				    EG_localToGlobal(tess, face+1, node+1,
						     &(nodes[0])), "l2g0");
			      nodes[0] -= 1;
			      printf("tess node index %d\n",nodes[0]);
			      ref_node_location(ref_grid_node(ref_grid),
						nodes[0]);
			      printf("u tess %f tmin %f tmax %f\n",
				     uv[0+2*node],uvmin[0],uvmax[0]);
			      printf("v tess %f tmin %f tmax %f\n",
				     uv[1+2*node],uvmin[1],uvmax[1]);
			      param[0] = 0.5*(uvmin[0]+uvmax[0]);
			      param[1] = 0.5*(uvmin[1]+uvmax[1]);
      RSS( ref_geom_add( ref_geom, nodes[0], REF_GEOM_FACE, face+1, param),
	   "face uv");

			    }
			}
		    }
		}
	    }
	}
    }

  }
  
  for (edge = 0; edge < (ref_geom->nedge); edge++) {
    int egads_status;
    REF_BOOL degenerate;
    degenerate = REF_FALSE;
    REIS( EGADS_SUCCESS,
	  EG_getTessEdge(tess, edge+1, &plen, &points, &t), "tess query edge" );
    for ( node = 0; node<plen; node++ ) {
      egads_status = EG_localToGlobal(tess, -(edge+1), node+1, &(nodes[0]));
      if ( EGADS_DEGEN == egads_status ) {
	degenerate = REF_TRUE;
      }else{
	REIS( EGADS_SUCCESS, egads_status, "l2g0");
	nodes[0] -= 1;
	param[0] = t[node];
	RSS( ref_geom_add( ref_geom, nodes[0], REF_GEOM_EDGE, edge+1, param),
	     "edge t");
      }
    }
    if ( !degenerate )
      for ( node = 0; node<(plen-1); node++ ) {
	/* assue edge index is 1-bias */
	REIS( EGADS_SUCCESS,
	      EG_localToGlobal(tess, -(edge+1), node+1, &(nodes[0])), "l2g0");
	REIS( EGADS_SUCCESS,
	      EG_localToGlobal(tess, -(edge+1), node+2, &(nodes[1])), "l2g1");
	nodes[0] -= 1;
	nodes[1] -= 1;
	nodes[2] = edge + 1;
	RSS( ref_cell_add(ref_grid_edg(ref_grid), nodes, &new_cell ),
	     "new edge");
      }
  }
  
#else
  printf("returning empty grid from %s, No EGADS linked.\n",__func__);
  SUPRESS_UNUSED_COMPILER_WARNING(ref_grid);
  SUPRESS_UNUSED_COMPILER_WARNING(max_length);
#endif
  
  return REF_SUCCESS;
}

REF_STATUS ref_geom_edge_tec_zone( REF_GRID ref_grid, REF_INT id, FILE *file )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_edg(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, nedg;

  RSS( ref_dict_create( &ref_dict ), "create dict" ); 

  each_ref_geom_edge( ref_geom, geom )
    if ( id == ref_geom_id(ref_geom,geom) )
      RSS( ref_dict_store( ref_dict, ref_geom_node(ref_geom,geom), geom ),
	   "mark nodes");
  nnode = ref_dict_n( ref_dict );
  
  nedg = 0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    if ( id == nodes[2] )
      nedg++;

  if ( 0 == nnode || 0 == nedg ) /* skip degenerate */
    {
      RSS( ref_dict_free( ref_dict ), "free dict" ); 
      return REF_SUCCESS;
    }
  
  fprintf(file,
	  "zone t=edge%d, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  id, nnode, nedg, "point", "felineseg" );

  each_ref_dict_key_value( ref_dict, item, node, geom )
    {
      fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
	      ref_node_xyz(ref_node,0,node),
	      ref_node_xyz(ref_node,1,node),
	      ref_node_xyz(ref_node,2,node),
	      0.0,
	      0.0,
	      ref_geom_param(ref_geom,0,geom) ) ;
    }

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( id == nodes[2] )
      {
	for ( node = 0; node < 2; node++ )
	  {
	    RSS( ref_dict_location( ref_dict, nodes[node], &local),
		 "localize");
	    fprintf(file," %d",local+1);
	  }
	fprintf(file,"\n");
      }

  RSS( ref_dict_free( ref_dict ), "free dict" ); 

  return REF_SUCCESS;
}

REF_STATUS ref_geom_face_tec_zone( REF_GRID ref_grid, REF_INT id, FILE *file )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, ntri;

  RSS( ref_dict_create( &ref_dict ), "create dict" ); 

  each_ref_geom_face( ref_geom, geom )
    if ( id == ref_geom_id(ref_geom,geom) )
      RSS( ref_dict_store( ref_dict, ref_geom_node(ref_geom,geom), geom ),
	   "mark nodes");
  nnode = ref_dict_n( ref_dict );
  
  ntri = 0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    if ( id == nodes[3] )
      ntri++;

  if ( 0 == nnode || 0 == ntri )  /* skip degenerate */
    {
      RSS( ref_dict_free( ref_dict ), "free dict" ); 
      return REF_SUCCESS;
    }
  
  fprintf(file,
	  "zone t=face%d, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  id, nnode, ntri, "point", "fetriangle" );

  each_ref_dict_key_value( ref_dict, item, node, geom )
    {
      fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
	      ref_node_xyz(ref_node,0,node),
	      ref_node_xyz(ref_node,1,node),
	      ref_node_xyz(ref_node,2,node),
	      ref_geom_param(ref_geom,0,geom),
	      ref_geom_param(ref_geom,1,geom),
	      0.0 ) ;
    }

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( id == nodes[3] )
      {
	for ( node = 0; node < 3; node++ )
	  {
	    RSS( ref_dict_location( ref_dict, nodes[node], &local),
		 "localize");
	    fprintf(file," %d",local+1);
	  }
	fprintf(file,"\n");
      }

  RSS( ref_dict_free( ref_dict ), "free dict" ); 

  return REF_SUCCESS;
}

REF_STATUS ref_geom_norm_tec_zone( REF_GRID ref_grid, REF_INT id, FILE *file )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, ntri;
  REF_DBL r[3], s[3], n[3];
  REF_DBL area_sign;

  RSS( ref_dict_create( &ref_dict ), "create dict" ); 

  each_ref_geom_face( ref_geom, geom )
    if ( id == ref_geom_id(ref_geom,geom) )
      RSS( ref_dict_store( ref_dict, ref_geom_node(ref_geom,geom), geom ),
	   "mark nodes");
  nnode = ref_dict_n( ref_dict );
  
  ntri = 0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    if ( id == nodes[3] )
      ntri++;

  if ( 0 == nnode || 0 == ntri ) /* skip degenerate */
    {
      RSS( ref_dict_free( ref_dict ), "free dict" ); 
      return REF_SUCCESS;
    }
  
  fprintf(file,
	  "zone t=norm%d, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  id, nnode, ntri, "point", "fetriangle" );
  
  each_ref_dict_key_value( ref_dict, item, node, geom )
    {
      RSS( ref_geom_find( ref_geom, node, REF_GEOM_FACE, id, &geom ),
	   "not found");
      RSS( ref_geom_rsn( ref_geom, geom, r, s, n ), "rsn");
      RSS( ref_geom_uv_area_sign( ref_grid, id, &area_sign ), "a sign" );
      fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
	      ref_node_xyz(ref_node,0,node),
	      ref_node_xyz(ref_node,1,node),
	      ref_node_xyz(ref_node,2,node),
	      -area_sign*n[0],
	      -area_sign*n[1],
	      -area_sign*n[2] ) ;
    }

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( id == nodes[3] )
      {
	for ( node = 0; node < 3; node++ )
	  {
	    RSS( ref_dict_location( ref_dict, nodes[node], &local),
		 "localize");
	    fprintf(file," %d",local+1);
	  }
	fprintf(file,"\n");
      }

  RSS( ref_dict_free( ref_dict ), "free dict" ); 

  return REF_SUCCESS;
}

REF_STATUS ref_geom_curve_tec_zone( REF_GRID ref_grid, REF_INT id, FILE *file )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_DICT ref_dict;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, ntri;
  REF_DBL kr, r[3], ks, s[3];

  RSS( ref_dict_create( &ref_dict ), "create dict" ); 

  each_ref_geom_face( ref_geom, geom )
    if ( id == ref_geom_id(ref_geom,geom) )
      RSS( ref_dict_store( ref_dict, ref_geom_node(ref_geom,geom), geom ),
	   "mark nodes");
  nnode = ref_dict_n( ref_dict );
  
  ntri = 0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    if ( id == nodes[3] )
      ntri++;

  if ( 0 == nnode || 0 == ntri ) /* skip degenerate */
    {
      RSS( ref_dict_free( ref_dict ), "free dict" ); 
      return REF_SUCCESS;
    }
  
  fprintf(file,
	  "zone t=curve%d, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  id, nnode, ntri, "point", "fetriangle" );
  
  each_ref_dict_key_value( ref_dict, item, node, geom )
    {
      RSS( ref_geom_find( ref_geom, node, REF_GEOM_FACE, id, &geom ),
	   "not found");
      RSS( ref_geom_curvature( ref_geom, geom, &kr, r, &ks, s ), "curve");
      fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
	      ref_node_xyz(ref_node,0,node),
	      ref_node_xyz(ref_node,1,node),
	      ref_node_xyz(ref_node,2,node),
	      ABS(kr),
	      ABS(ks),
	      0.0) ;
    }

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( id == nodes[3] )
      {
	for ( node = 0; node < 3; node++ )
	  {
	    RSS( ref_dict_location( ref_dict, nodes[node], &local),
		 "localize");
	    fprintf(file," %d",local+1);
	  }
	fprintf(file,"\n");
      }

  RSS( ref_dict_free( ref_dict ), "free dict" ); 

  return REF_SUCCESS;
}

REF_STATUS ref_geom_tec( REF_GRID ref_grid, const char *filename  )
{
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  FILE *file;
  REF_INT geom, id, min_id, max_id;
  
  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  fprintf(file, "title=\"refine cad coupling in tecplot format\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\" \"u\" \"v\" \"t\"\n");

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_edge( ref_geom, geom )
    {
      min_id = MIN( min_id, ref_geom_id(ref_geom,geom) );
      max_id = MAX( max_id, ref_geom_id(ref_geom,geom) );
    }

  for ( id = min_id ; id <= max_id ; id++ )
    RSS( ref_geom_edge_tec_zone( ref_grid, id, file ), "tec edge" );

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_face( ref_geom, geom )
    {
      min_id = MIN( min_id, ref_geom_id(ref_geom,geom) );
      max_id = MAX( max_id, ref_geom_id(ref_geom,geom) );
    }
  
  for ( id = min_id ; id <= max_id ; id++ )
    RSS( ref_geom_face_tec_zone( ref_grid, id, file ), "tec face" );
  for ( id = min_id ; id <= max_id ; id++ )
    RSS( ref_geom_norm_tec_zone( ref_grid, id, file ), "tec norm" );
  for ( id = min_id ; id <= max_id ; id++ )
    RSS( ref_geom_curve_tec_zone( ref_grid, id, file ), "tec curve" );

  fclose(file);
  return REF_SUCCESS;
}

REF_STATUS ref_geom_ghost( REF_GEOM ref_geom, REF_NODE ref_node )
{
  REF_INT *a_nnode, *b_nnode;
  REF_INT a_nnode_total, b_nnode_total;
  REF_INT *a_global, *b_global;
  REF_INT *a_part, *b_part;
  REF_INT *a_ngeom, *b_ngeom;
  REF_INT a_ngeom_total, b_ngeom_total;
  REF_INT *a_tgi, *b_tgi;
  REF_DBL *a_param, *b_param;
  REF_INT part, node, degree;
  REF_INT *a_next, *b_next;
  REF_INT local, item, geom;

  if ( 1 == ref_mpi_n ) return REF_SUCCESS;

  ref_malloc_init( a_next,  ref_mpi_n, REF_INT, 0 );
  ref_malloc_init( b_next,  ref_mpi_n, REF_INT, 0 );
  ref_malloc_init( a_nnode, ref_mpi_n, REF_INT, 0 );
  ref_malloc_init( b_nnode, ref_mpi_n, REF_INT, 0 );
  ref_malloc_init( a_ngeom, ref_mpi_n, REF_INT, 0 );
  ref_malloc_init( b_ngeom, ref_mpi_n, REF_INT, 0 );

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      a_nnode[ref_node_part(ref_node,node)]++;

  RSS( ref_mpi_alltoall( a_nnode, b_nnode, REF_INT_TYPE ), "alltoall nnodes");

  a_nnode_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    a_nnode_total += a_nnode[part];
  ref_malloc( a_global, a_nnode_total, REF_INT );
  ref_malloc( a_part, a_nnode_total, REF_INT );

  b_nnode_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    b_nnode_total += b_nnode[part];
  ref_malloc( b_global, b_nnode_total, REF_INT );
  ref_malloc( b_part, b_nnode_total, REF_INT );

  a_next[0] = 0;
  for ( part = 1; part<ref_mpi_n ; part++ )
    a_next[part] = a_next[part-1]+a_nnode[part-1];

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      {
	part = ref_node_part(ref_node,node);
	a_global[a_next[part]] = ref_node_global(ref_node,node);
	a_part[a_next[part]] = ref_mpi_id;
	a_next[ref_node_part(ref_node,node)]++;
      }

  RSS( ref_mpi_alltoallv( a_global, a_nnode, b_global, b_nnode, 
			  1, REF_INT_TYPE ), 
       "alltoallv global");
  RSS( ref_mpi_alltoallv( a_part, a_nnode, b_part, b_nnode, 
			  1, REF_INT_TYPE ), 
       "alltoallv global");

  for (node=0;node<b_nnode_total;node++)
    {
      RSS( ref_node_local( ref_node, b_global[node], &local ), "g2l");
      part = b_part[node];
      RSS( ref_adj_degree( ref_geom_adj(ref_geom), local, &degree ), "deg" );
      /* printf("%d: node %d global %d local %d part %d degree %d\n",
	 ref_mpi_id, node,b_global[node], local, part, degree); */
      b_ngeom[part] += degree;
    }

  RSS( ref_mpi_alltoall( b_ngeom, a_ngeom, REF_INT_TYPE ), "alltoall ngeoms");

  a_ngeom_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    a_ngeom_total += a_ngeom[part];
  ref_malloc( a_tgi,   3*a_ngeom_total, REF_INT );
  ref_malloc( a_param, 2*a_ngeom_total, REF_DBL );

  b_ngeom_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    b_ngeom_total += b_ngeom[part];
  ref_malloc( b_tgi,   3*b_ngeom_total, REF_INT );
  ref_malloc( b_param, 2*b_ngeom_total, REF_DBL );

  b_next[0] = 0;
  for ( part = 1; part<ref_mpi_n ; part++ )
    b_next[part] = b_next[part-1]+b_ngeom[part-1];

  for (node=0;node<b_nnode_total;node++)
    {
      RSS( ref_node_local( ref_node, b_global[node], &local ), "g2l");
      part = b_part[node];
      each_ref_geom_having_node( ref_geom, local, item, geom )
	{
	  b_tgi[0+3*b_next[part]] = ref_geom_type(ref_geom,geom);
	  b_tgi[1+3*b_next[part]] =
	    ref_node_global(ref_node,ref_geom_node(ref_geom,geom));
	  b_tgi[2+3*b_next[part]] = ref_geom_id(ref_geom,geom);
	  b_param[0+2*b_next[part]] = ref_geom_param(ref_geom,0,geom);
	  b_param[1+2*b_next[part]] = ref_geom_param(ref_geom,1,geom);
	  b_next[part]++;
	}
    }

  RSS( ref_mpi_alltoallv( b_tgi, b_ngeom, a_tgi, a_ngeom, 
			  3, REF_INT_TYPE ), 
       "alltoallv tgi");
  RSS( ref_mpi_alltoallv( b_param, b_ngeom, a_param, a_ngeom, 
			  2, REF_DBL_TYPE ), 
       "alltoallv param");

  for (geom=0;geom<a_ngeom_total;geom++)
    {
      RSS( ref_node_local( ref_node, a_tgi[1+3*geom], &local ), "g2l");
      RSS( ref_geom_add( ref_geom, local, a_tgi[0+3*geom], a_tgi[2+3*geom],
			 &(a_param[2*geom]) ), "add ghost" );
    }
  
  free(b_param);
  free(b_tgi);
  free(a_param);
  free(a_tgi);
  free(b_part);
  free(b_global);
  free(a_part);
  free(a_global);
  free(b_ngeom);
  free(a_ngeom);
  free(b_nnode);
  free(a_nnode);
  free(b_next);
  free(a_next);

  return REF_SUCCESS;  
}

REF_STATUS ref_geom_faceid_range( REF_GRID ref_grid, 
				  REF_INT *min_faceid, REF_INT *max_faceid )
{
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  *min_faceid = REF_INT_MAX;
  *max_faceid = REF_INT_MIN;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      *min_faceid = MIN( *min_faceid, nodes[ref_cell_node_per(ref_cell)] );
      *max_faceid = MAX( *max_faceid, nodes[ref_cell_node_per(ref_cell)] );
    }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      *min_faceid = MIN( *min_faceid, nodes[ref_cell_node_per(ref_cell)] );
      *max_faceid = MAX( *max_faceid, nodes[ref_cell_node_per(ref_cell)] );
    }

  if ( ref_mpi_n > 1 )
    {
      REF_INT global;

      RSS( ref_mpi_min( min_faceid, &global, REF_INT_TYPE ), "mpi min face" );
      RSS( ref_mpi_bcast( &global, 1, REF_INT_TYPE ), "mpi min face" );
      *min_faceid = global;

      RSS( ref_mpi_max( max_faceid, &global, REF_INT_TYPE ), "mpi max face" );
      RSS( ref_mpi_bcast( &global, 1, REF_INT_TYPE ), "mpi max face" );
      *max_faceid = global;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_geom_edgeid_range( REF_GRID ref_grid, 
				  REF_INT *min_edgeid, REF_INT *max_edgeid )
{
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  *min_edgeid = REF_INT_MAX;
  *max_edgeid = REF_INT_MIN;

  ref_cell = ref_grid_edg(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      *min_edgeid = MIN( *min_edgeid, nodes[ref_cell_node_per(ref_cell)] );
      *max_edgeid = MAX( *max_edgeid, nodes[ref_cell_node_per(ref_cell)] );
    }

  if ( ref_mpi_n > 1 )
    {
      REF_INT global;

      RSS( ref_mpi_min( min_edgeid, &global, REF_INT_TYPE ), "mpi min edge" );
      RSS( ref_mpi_bcast( &global, 1, REF_INT_TYPE ), "mpi min edge" );
      *min_edgeid = global;

      RSS( ref_mpi_max( max_edgeid, &global, REF_INT_TYPE ), "mpi max edge" );
      RSS( ref_mpi_bcast( &global, 1, REF_INT_TYPE ), "mpi max edge" );
      *max_edgeid = global;
    }

  return REF_SUCCESS;
}

