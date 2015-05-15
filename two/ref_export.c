
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_export.h"

#include "ref_dict.h"
#include "ref_endian.h"

#include "ref_mpi.h"

#include "ref_matrix.h"

#include "ref_edge.h"
#include "ref_node.h"

#include "ref_malloc.h"
#include "ref_math.h"

#define VTK_TETRA      (10)
#define VTK_HEXAHEDRON (12)
#define VTK_WEDGE      (13)
#define VTK_PYRAMID    (14)

/*
3-4 UGRID
| |\ 
| | 2
| |/
0-1
2-3 VTK
| |\ 
| | 4
| |/
1-0
 */

#define VTK_PYRAMID_ORDER(vtk_nodes)			\
  {							\
    REF_INT ugrid_nodes[5];				\
    ugrid_nodes[0] = (vtk_nodes)[0];			\
    ugrid_nodes[1] = (vtk_nodes)[1];			\
    ugrid_nodes[2] = (vtk_nodes)[2];			\
    ugrid_nodes[3] = (vtk_nodes)[3];			\
    ugrid_nodes[4] = (vtk_nodes)[4];			\
    (vtk_nodes)[0] = ugrid_nodes[1];			\
    (vtk_nodes)[1] = ugrid_nodes[0];			\
    (vtk_nodes)[2] = ugrid_nodes[3];			\
    (vtk_nodes)[3] = ugrid_nodes[4];			\
    (vtk_nodes)[4] = ugrid_nodes[2];			\
  }

/*
  tecplot "brick" 
      7---6
     /|  /|
    4-+-5 |
    | | | |
    | 3-+-2
    |/  |/
    0---1
 */

#define TEC_BRICK_TET(brick,nodes)					\
  {									\
    brick[0] = nodes[0]; brick[1] = nodes[1]; brick[2] = nodes[2];	\
    brick[3] = nodes[2];						\
    brick[4] = nodes[3]; brick[5] = nodes[3];				\
    brick[6] = nodes[3]; brick[7] = nodes[3];				\
  }

#define TEC_BRICK_PYR(brick,nodes)					\
  {									\
    brick[0] = nodes[0]; brick[1] = nodes[1];				\
    brick[2] = nodes[2]; brick[3] = nodes[2];				\
    brick[4] = nodes[3]; brick[5] = nodes[4];				\
    brick[6] = nodes[2]; brick[7] = nodes[2];				\
  }

#define TEC_BRICK_PRI(brick,nodes)					\
  {									\
    brick[0] = nodes[0]; brick[1] = nodes[1];				\
    brick[2] = nodes[2]; brick[3] = nodes[2];				\
    brick[4] = nodes[3]; brick[5] = nodes[4];				\
    brick[6] = nodes[5]; brick[7] = nodes[5];				\
  }

#define TEC_BRICK_HEX(brick,nodes)					\
  {									\
    brick[0] = nodes[0]; brick[1] = nodes[1];				\
    brick[2] = nodes[2]; brick[3] = nodes[3];				\
    brick[4] = nodes[4]; brick[5] = nodes[5];				\
    brick[6] = nodes[6]; brick[7] = nodes[7];				\
  }

REF_STATUS ref_export_by_extension( REF_GRID ref_grid, char *filename )
{
  size_t end_of_string;

  end_of_string = strlen(filename);

  if( strcmp(&filename[end_of_string-4],".vtk") == 0 ) 
    {
      RSS( ref_export_vtk( ref_grid, filename ), "vtk export failed");
    } 
  else 
    if( strcmp(&filename[end_of_string-2],".c") == 0 ) 
      {
	RSS( ref_export_c( ref_grid, filename ), "C export failed");
      } 
    else 
      if( strcmp(&filename[end_of_string-4],".tec") == 0 ) 
	{
	  RSS( ref_export_tec( ref_grid, filename ), "tec export failed");
	} 
      else 
	if( strcmp(&filename[end_of_string-4],".eps") == 0 ) 
	  {
	    RSS( ref_export_eps( ref_grid, filename ), "eps export failed");
	  } 
	else 
	if( strcmp(&filename[end_of_string-4],".pdf") == 0 ) 
	  {
	    RSS( ref_export_pdf( ref_grid, filename ), "pdf export failed");
	  } 
	else 
	  if( strcmp(&filename[end_of_string-10],".lb8.ugrid") == 0 ) 
	    {
	      RSS( ref_export_lb8_ugrid( ref_grid, filename ), 
		   "lb8.ugrid export failed");
	    } 
	  else 
	  if( strcmp(&filename[end_of_string-9],".b8.ugrid") == 0 ) 
	    {
	      RSS( ref_export_b8_ugrid( ref_grid, filename ), 
		   "b8.ugrid export failed");
	    } 
	  else 
	    if( strcmp(&filename[end_of_string-6],".ugrid") == 0 ) 
	      {
		RSS( ref_export_ugrid( ref_grid, filename ), 
		     "ugrid export failed");
	      } 
	    else 
	      if( strcmp(&filename[end_of_string-6],".fgrid") == 0 ) 
		{
		  RSS( ref_export_fgrid( ref_grid, filename ), 
		       "ugrid export failed");
		} 
	      else 
		if( strcmp(&filename[end_of_string-6],".cogsg") == 0 ) 
		  {
		    RSS( ref_export_cogsg( ref_grid, filename ), 
			 "cogsg export failed");
		  } 
		else 
		  if( strcmp(&filename[end_of_string-5],".html") == 0 ) 
		    {
		      RSS( ref_export_html( ref_grid, filename ), 
			   "html export failed");
		    } 
		  else 
		    if( strcmp(&filename[end_of_string-9],".2d.meshb") == 0 ) 
		      {
			RSS( ref_export_twod_meshb( ref_grid, filename ), 
			     "twod meshb export failed");
		      } 
		  else 
		    if( strcmp(&filename[end_of_string-6],".meshb") == 0 ) 
		      {
			RSS( ref_export_meshb( ref_grid, filename ), 
			     "meshb export failed");
		      } 
		  else 
		    if( strcmp(&filename[end_of_string-4],".msh") == 0 ) 
		      {
			RSS( ref_export_twod_msh( ref_grid, filename ), 
			     "msh export failed");
		      } 
		    else 
		      {
			printf("%s: %d: %s %s\n",__FILE__,__LINE__,
			       "export file name extension unknown", filename);
			RSS( REF_FAILURE, "unknown file extension");
		      }

  return REF_SUCCESS;
}

REF_STATUS ref_export_vtk( REF_GRID ref_grid, char *filename  )
{
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT ncell,size;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell;
  REF_INT group;

  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  fprintf(file,"# vtk DataFile Version 2.0\n");
  fprintf(file,"ref_export_vtk\n");
  fprintf(file,"ASCII\n");

  RSS( ref_node_compact( ref_node, &o2n, &n2o ), "compact" );

  fprintf(file,"DARASET UNSTRUCTURED_GRID\n");
  fprintf(file,"POINTS %d double\n",ref_node_n(ref_node));

  for ( node = 0; node < ref_node_n(ref_node); node++ )
    {
      fprintf(file, " %.16e %.16e %.16e\n",
	      ref_node_xyz(ref_node,0,n2o[node]),
	      ref_node_xyz(ref_node,1,n2o[node]),
	      ref_node_xyz(ref_node,2,n2o[node]) ) ;
    }

  ncell = 0;
  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    ncell += ref_cell_n(ref_cell);

  size = 0;
  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    size += ref_cell_n(ref_cell)*(1+ref_cell_node_per(ref_cell));

  fprintf(file,"CELLS %d %d\n",ncell,size);

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    {
      node_per = ref_cell_node_per(ref_cell);
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  fprintf(file," %d",node_per);
	  if ( 5 == node_per ) VTK_PYRAMID_ORDER(nodes);
	  for ( node = 0; node < node_per; node++ )
	    fprintf(file," %d",o2n[nodes[node]]);
	  fprintf(file,"\n");
	}
    }

  fprintf(file,"CELL_TYPES %d\n",ncell);

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_valid_cell( ref_cell, cell )
    fprintf(file," %d\n",VTK_TETRA);

  ref_cell = ref_grid_pyr(ref_grid);
  each_ref_cell_valid_cell( ref_cell, cell )
    fprintf(file," %d\n",VTK_PYRAMID);

  ref_cell = ref_grid_pri(ref_grid);
  each_ref_cell_valid_cell( ref_cell, cell )
    fprintf(file," %d\n",VTK_WEDGE);

  ref_cell = ref_grid_hex(ref_grid);
  each_ref_cell_valid_cell( ref_cell, cell )
    fprintf(file," %d\n",VTK_HEXAHEDRON);

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec( REF_GRID ref_grid, char *filename  )
{
  FILE *file;

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  fprintf(file, "title=\"tecplot refine geometry file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

  RSS( ref_export_tec_surf_zone( ref_grid, file  ), "surf" );
  RSS( ref_export_tec_vol_zone( ref_grid, file  ), "vol" );

  fclose(file);
  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_surf( REF_GRID ref_grid, char *filename  )
{
  FILE *file;

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  fprintf(file, "title=\"tecplot refine geometry file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

  RSS( ref_export_tec_surf_zone( ref_grid, file  ), "surf" );

  fclose(file);
  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_surf_zone( REF_GRID ref_grid, FILE *file )
{
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *g2l, *l2g;
  REF_INT nface;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnode;
  REF_DICT ref_dict;
  REF_INT boundary_tag,boundary_index;

  ref_node = ref_grid_node(ref_grid);

  RSS( ref_dict_create( &ref_dict ), "create dict" ); 

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    RSS( ref_dict_store( ref_dict, nodes[3], REF_EMPTY ), "mark tri" );

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    RSS( ref_dict_store( ref_dict, nodes[4], REF_EMPTY ), "mark qua" );

  each_ref_dict_key( ref_dict, boundary_index, boundary_tag )
    {
      RSS( ref_grid_boundary_nodes( ref_grid, boundary_tag,
				    &nnode, &nface, &g2l, &l2g ),
	   "extract this boundary");

      fprintf(file,
	  "zone t=surf%d, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	      boundary_tag, nnode, nface, "point", "fequadrilateral" );

      for ( node = 0; node < nnode; node++ )
	fprintf(file, " %.16e %.16e %.16e\n",
		ref_node_xyz(ref_node,0,l2g[node]),
		ref_node_xyz(ref_node,1,l2g[node]),
		ref_node_xyz(ref_node,2,l2g[node]) ) ;

      ref_cell = ref_grid_tri(ref_grid);
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	if ( boundary_tag == nodes[3] )
	  {
	    nodes[3] = nodes[2];
	    for ( node = 0; node < 4; node++ )
	      {
		fprintf(file," %d",g2l[nodes[node]] + 1);
	      }
	    fprintf(file,"\n");
	  }

      ref_cell = ref_grid_qua(ref_grid);
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	if ( boundary_tag == nodes[4] )
	  {
	    for ( node = 0; node < 4; node++ )
	      fprintf(file," %d",g2l[nodes[node]] + 1);
	    fprintf(file,"\n");
	  }

      ref_free(l2g);
      ref_free(g2l);
    }

  RSS( ref_dict_free( ref_dict ), "free dict" ); 

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_vol_zone( REF_GRID ref_grid, FILE *file  )
{
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT brick[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnode;
  REF_INT group, node_per;

  ref_node = ref_grid_node(ref_grid);

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    if ( ref_cell_n(ref_cell) > 0 )
      {
	node_per = ref_cell_node_per(ref_cell);

	o2n = (REF_INT *)malloc( ref_node_max(ref_node) * sizeof(REF_INT) );
	RNS(o2n,"malloc o2n NULL");

	nnode = 0;
	for ( node = 0 ; node < ref_node_max(ref_node) ; node++ )
	  o2n[node] = REF_EMPTY;

	each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	  for ( node = 0; node < node_per; node++ )
	    if ( REF_EMPTY == o2n[nodes[node]] )
	      { o2n[nodes[node]] = nnode; nnode++; }

	n2o = (REF_INT *)malloc( nnode * sizeof(REF_INT) );
	RNS(n2o,"malloc n2o NULL");

	for ( node = 0 ; node < ref_node_max(ref_node) ; node++ )
	  if ( REF_EMPTY != o2n[node] ) n2o[o2n[node]] = node;

	fprintf(file,
		"zone t=e%d, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
		node_per, nnode, ref_cell_n(ref_cell), "point", "febrick" );

	for ( node = 0; node < nnode; node++ )
	  fprintf(file, " %.16e %.16e %.16e\n",
		  ref_node_xyz(ref_node,0,n2o[node]),
		  ref_node_xyz(ref_node,1,n2o[node]),
		  ref_node_xyz(ref_node,2,n2o[node]) ) ;

	each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	  {
	    switch ( ref_cell_node_per(ref_cell) )
	      {
	      case 4:
		TEC_BRICK_TET(brick,nodes);
		break;
	      case 5:
		TEC_BRICK_PYR(brick,nodes);
		break;
	      case 6:
		TEC_BRICK_PRI(brick,nodes);
		break;
	      case 8:
		TEC_BRICK_HEX(brick,nodes);
		break;
	      default:
		RSS( REF_IMPLEMENT, "wrong nodes per cell");
		break;
	      }

	    for ( node = 0; node < 8; node++ )
	      {
		fprintf(file," %d",o2n[brick[node]] + 1);
	      }
	    fprintf(file,"\n");
	  }

	ref_free(n2o);
	ref_free(o2n);

      }

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_int( REF_GRID ref_grid, REF_INT *scalar, 
			       char *filename  )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT brick[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT ncell;
  REF_INT group;

  FILE *file;

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  fprintf(file, "title=\"tecplot refine scalar file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\" \"s\"\n");

  ncell = 0;
  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    ncell += ref_cell_n(ref_cell);

  fprintf(file,
	  "zone t=scalar, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  ref_node_n(ref_node), ncell, "point", "febrick" );

  RSS( ref_node_compact( ref_node, &o2n, &n2o), "compact" );

  for ( node = 0; node < ref_node_n(ref_node); node++ )
    fprintf(file, " %.16e %.16e %.16e %d\n", 
	    ref_node_xyz(ref_node,0,n2o[node]),
	    ref_node_xyz(ref_node,1,n2o[node]),
	    ref_node_xyz(ref_node,2,n2o[node]),
	    scalar[n2o[node]] ) ;

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      switch ( ref_cell_node_per(ref_cell) )
	{
	case 4:
	  TEC_BRICK_TET(brick,nodes);
	  break;
	case 5:
	  TEC_BRICK_PYR(brick,nodes);
	  break;
	case 6:
	  TEC_BRICK_PRI(brick,nodes);
	  break;
	case 8:
	  TEC_BRICK_HEX(brick,nodes);
	  break;
	default:
	  RSS( REF_IMPLEMENT, "wrong nodes per cell");
	  break;
	}

      for ( node = 0; node < 8; node++ )
	{
	  fprintf(file," %d",o2n[brick[node]] + 1);
	}
      fprintf(file,"\n");
    }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_part( REF_GRID ref_grid, char *root_filename )
{
  REF_NODE ref_node = ref_grid_node( ref_grid );
  char viz_file[256];

  sprintf(viz_file, "%s_n%d_p%d.tec", root_filename, ref_mpi_n, ref_mpi_id);

  RSS(ref_export_tec_int( ref_grid, ref_node->part,
			  viz_file ) , "viz parts as scalar");

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_metric_axis( REF_GRID ref_grid, char *root_filename )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT brick[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT ncell;
  REF_INT group;
  REF_DBL d[12];
  REF_DBL dx,dy,dz;
  FILE *file;
  char viz_file[256];
  REF_INT e;

  for (e=0;e<3;e++)
    {
      sprintf(viz_file, "%s_n%d_p%d_ax%d.tec", 
	      root_filename, ref_mpi_n, ref_mpi_id, e);

      file = fopen(viz_file,"w");
      if (NULL == (void *)file) printf("unable to open %s\n",viz_file);
      RNS(file, "unable to open file" );

      fprintf(file, "title=\"tecplot refine metric axes\"\n");
      fprintf(file, "variables = \"x\" \"y\" \"z\" \"u\" \"v\" \"w\"\n");

      ncell = 0;
      each_ref_grid_ref_cell( ref_grid, group, ref_cell )
	ncell += ref_cell_n(ref_cell);

      fprintf(file,
	  "zone t=scalar, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  ref_node_n(ref_node), ncell, "point", "febrick" );

      RSS( ref_node_compact( ref_node, &o2n, &n2o ), "compact" );

      for ( node = 0; node < ref_node_n(ref_node); node++ )
	{
	  RSS( ref_matrix_diag_m( ref_node_metric_ptr(ref_node,n2o[node]),
				  d ), "diag" );
	  RSS( ref_matrix_ascending_eig( d ), "sort eig" );
	  dx = d[3+3*e]/sqrt(d[e]);
	  dy = d[4+3*e]/sqrt(d[e]);
	  dz = d[5+3*e]/sqrt(d[e]);
	  fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e\n", 
		  ref_node_xyz(ref_node,0,n2o[node]),
		  ref_node_xyz(ref_node,1,n2o[node]),
		  ref_node_xyz(ref_node,2,n2o[node]),
		  dx, dy, dz);
	}

      each_ref_grid_ref_cell( ref_grid, group, ref_cell )
	each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  switch ( ref_cell_node_per(ref_cell) )
	    {
	    case 4:
	      TEC_BRICK_TET(brick,nodes);
	      break;
	    case 5:
	      TEC_BRICK_PYR(brick,nodes);
	      break;
	    case 6:
	      TEC_BRICK_PRI(brick,nodes);
	      break;
	    case 8:
	      TEC_BRICK_HEX(brick,nodes);
	      break;
	    default:
	      RSS( REF_IMPLEMENT, "wrong nodes per cell");
	      break;
	    }

	  for ( node = 0; node < 8; node++ )
	    {
	      fprintf(file," %d",o2n[brick[node]] + 1);
	    }
	  fprintf(file,"\n");
	}

      ref_free(n2o);
      ref_free(o2n);

      fclose(file);
    } /* each eigenpair */

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_metric_ellipse( REF_GRID ref_grid, 
					  char *root_filename )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT ncell;
  REF_DBL d[12];
  REF_DBL x,y,z;
  REF_DBL ex,ey;
  FILE *file;
  char viz_file[256];
  REF_INT i, n=36;
  REF_INT e0, e1;
  REF_DBL dt = ref_math_in_radians(360.0/(REF_DBL)n);
  REF_DBL scale = 0.5; /* so the ellipses touch for an ideal grid */

  sprintf(viz_file, "%s_n%d_p%d_ellipse.tec", 
	  root_filename, ref_mpi_n, ref_mpi_id);

  file = fopen(viz_file,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",viz_file);
  RNS(file, "unable to open file" );

  fprintf(file, "title=\"tecplot refine metric axes\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

  ncell = ref_node_n(ref_node)*n;

  fprintf(file,
	  "zone t=scalar, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  3*ncell, 3*ncell, "point", "felineseg" );

  RSS( ref_node_compact( ref_node, &o2n, &n2o ), "compact" );

  for ( node = 0; node < ref_node_n(ref_node); node++ )
    {
      RSS( ref_matrix_diag_m( ref_node_metric_ptr(ref_node,n2o[node]),
			      d ), "diag" );
      RSS( ref_matrix_ascending_eig( d ), "sort eig" );
      for (e0=0;e0<3;e0++)
	{
	  e1 = e0+1;
	  if (e1==3) 
	    e1 = 0;
	  for (i=0;i<n;i++)
	    {
	      ex = scale*cos(i*dt)/sqrt(d[e0]);
	      ey = scale*sin(i*dt)/sqrt(d[e1]);
	      x = d[3+3*e0]*ex + d[3+3*e1]*ey;
	      y = d[4+3*e0]*ex + d[4+3*e1]*ey;
	      z = d[5+3*e0]*ex + d[5+3*e1]*ey;
	      fprintf(file, " %.16e %.16e %.16e\n", 
		      ref_node_xyz(ref_node,0,n2o[node])+x,
		      ref_node_xyz(ref_node,1,n2o[node])+y,
		      ref_node_xyz(ref_node,2,n2o[node])+z);
	    }
	}
    }
  
  for (e0=0;e0<3;e0++)
    for ( node = 0; node < ref_node_n(ref_node); node++ )
      {
	for (i=0;i<n-1;i++)
	  fprintf(file," %d %d\n",i+node*n+1+ncell*e0,i+1+node*n+1+ncell*e0);
	fprintf(file," %d %d\n",n+node*n+ncell*e0,1+node*n+ncell*e0);
      }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_metric_box( REF_GRID ref_grid, 
				      char *root_filename,
				      REF_DBL *bounding_box)
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT node, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT nnode, *o2n, *n2o;
  REF_INT ncell;
  REF_DBL d[12];
  REF_DBL x,y,z;
  REF_DBL ex,ey;
  FILE *file;
  char viz_file[256];
  REF_INT i, n=36;
  REF_INT e0, e1;
  REF_DBL dt = ref_math_in_radians(360.0/(REF_DBL)n);
  REF_DBL scale = 0.5; /* so the ellipses touch for an ideal grid */

  sprintf(viz_file, "%s_n%d_p%d_ellipse.tec", 
	  root_filename, ref_mpi_n, ref_mpi_id);

  file = fopen(viz_file,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",viz_file);
  RNS(file, "unable to open file" );

  fprintf(file, "title=\"tecplot refine metric axes\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

  RSS( ref_node_in_bounding_box( ref_node, bounding_box,
				 &nnode, &o2n, &n2o ), "bbox" );

  ncell = nnode*n;

  fprintf(file,
	  "zone t=met, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  3*ncell, 3*ncell, "point", "felineseg" );

  for ( node = 0; node < nnode; node++ )
    {
      RSS( ref_matrix_diag_m( ref_node_metric_ptr(ref_node,n2o[node]),
			      d ), "diag" );
      RSS( ref_matrix_ascending_eig( d ), "sort eig" );
      for (e0=0;e0<3;e0++)
	{
	  e1 = e0+1;
	  if (e1==3) 
	    e1 = 0;
	  for (i=0;i<n;i++)
	    {
	      ex = scale*cos(i*dt)/sqrt(d[e0]);
	      ey = scale*sin(i*dt)/sqrt(d[e1]);
	      x = d[3+3*e0]*ex + d[3+3*e1]*ey;
	      y = d[4+3*e0]*ex + d[4+3*e1]*ey;
	      z = d[5+3*e0]*ex + d[5+3*e1]*ey;
	      fprintf(file, " %.16e %.16e %.16e\n", 
		      ref_node_xyz(ref_node,0,n2o[node])+x,
		      ref_node_xyz(ref_node,1,n2o[node])+y,
		      ref_node_xyz(ref_node,2,n2o[node])+z);
	    }
	}
    }
  
  for (e0=0;e0<3;e0++)
    for ( node = 0; node < nnode; node++ )
      {
	for (i=0;i<n-1;i++)
	  fprintf(file," %d %d\n",i+node*n+1+ncell*e0,i+1+node*n+1+ncell*e0);
	fprintf(file," %d %d\n",n+node*n+ncell*e0,1+node*n+ncell*e0);
      }

  ncell = 0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( REF_EMPTY != o2n[nodes[0]] &&
	 REF_EMPTY != o2n[nodes[1]] &&
	 REF_EMPTY != o2n[nodes[2]] ) 
      ncell++;
 
  fprintf(file,
	  "zone t=tri, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  nnode, ncell, "point", "fetriangle" );

  for ( node = 0; node < nnode; node++ )
    {
      fprintf(file, " %.16e %.16e %.16e\n", 
	      ref_node_xyz(ref_node,0,n2o[node]),
	      ref_node_xyz(ref_node,1,n2o[node]),
	      ref_node_xyz(ref_node,2,n2o[node]));
    }

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( REF_EMPTY != o2n[nodes[0]] &&
	 REF_EMPTY != o2n[nodes[1]] &&
	 REF_EMPTY != o2n[nodes[2]] )
      {
	for ( node = 0; node < 3; node++ )
	  fprintf(file," %d",o2n[nodes[node]] + 1);
	fprintf(file,"\n");
      }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_ratio( REF_GRID ref_grid, char *root_filename )
{
  REF_NODE ref_node = ref_grid_node( ref_grid );
  REF_EDGE ref_edge;
  char viz_file[256];

  sprintf(viz_file, "%s_n%d_p%d.tec", root_filename, ref_mpi_n, ref_mpi_id);

  RSS( ref_edge_create( &ref_edge, ref_grid ), "make edge" );

  RSS(ref_edge_tec_ratio( ref_edge, ref_node,
			  viz_file ) , "viz parts as scalar");

  RSS( ref_edge_free( ref_edge ), "free edge" );

  return REF_SUCCESS;
}

REF_STATUS ref_export_fgrid( REF_GRID ref_grid, char *filename  )
{
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nnode,ntri,ntet;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell;
  REF_INT ixyz;

  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  nnode = ref_node_n(ref_node);

  ntri = ref_cell_n(ref_grid_tri(ref_grid));

  ntet = ref_cell_n(ref_grid_tet(ref_grid));

  fprintf(file,"%d %d %d\n",nnode,ntri,ntet);

  RSS( ref_node_compact( ref_node, &o2n, &n2o), "compact" );

  for ( ixyz = 0 ; ixyz< 3; ixyz++)
    for ( node = 0; node < ref_node_n(ref_node); node++ )
      fprintf(file, " %.16e\n", ref_node_xyz(ref_node,ixyz,n2o[node]) ) ;

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      for ( node = 0; node < node_per; node++ )
	fprintf(file," %d",o2n[nodes[node]]+1);
      fprintf(file,"\n");
    }
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      fprintf(file," %d",nodes[3]);
      fprintf(file,"\n");
    }

  ref_cell = ref_grid_tet(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      for ( node = 0; node < node_per; node++ )
	fprintf(file," %d",o2n[nodes[node]]+1);
      fprintf(file,"\n");
    }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_ugrid( REF_GRID ref_grid, char *filename  )
{
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nnode,ntri,nqua,ntet,npyr,npri,nhex;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell;
  REF_INT group;
  REF_INT faceid, min_faceid, max_faceid;

  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  nnode = ref_node_n(ref_node);

  ntri = ref_cell_n(ref_grid_tri(ref_grid));
  nqua = ref_cell_n(ref_grid_qua(ref_grid));

  ntet = ref_cell_n(ref_grid_tet(ref_grid));
  npyr = ref_cell_n(ref_grid_pyr(ref_grid));
  npri = ref_cell_n(ref_grid_pri(ref_grid));
  nhex = ref_cell_n(ref_grid_hex(ref_grid));

  fprintf(file,"%d %d %d %d %d %d %d\n",nnode,ntri,nqua,ntet,npyr,npri,nhex);

  RSS( ref_node_compact( ref_node, &o2n, &n2o), "compact" );

  for ( node = 0; node < ref_node_n(ref_node); node++ )
    fprintf(file, " %.16e %.16e %.16e\n", 
	    ref_node_xyz(ref_node,0,n2o[node]),
	    ref_node_xyz(ref_node,1,n2o[node]),
	    ref_node_xyz(ref_node,2,n2o[node]) ) ;

  RSS( ref_export_faceid_range( ref_grid, &min_faceid, &max_faceid), "range");

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      if ( nodes[node_per] == faceid )
	{
	  for ( node = 0; node < node_per; node++ )
	    fprintf(file," %d",o2n[nodes[node]]+1);
	  fprintf(file,"\n");
	}

  ref_cell = ref_grid_qua(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      if ( nodes[node_per] == faceid )
	{
	  for ( node = 0; node < node_per; node++ )
	    fprintf(file," %d",o2n[nodes[node]]+1);
	  fprintf(file,"\n");
	}

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      if ( nodes[node_per] == faceid )
	{
	  fprintf(file," %d",nodes[3]);
	  fprintf(file,"\n");
	}

  ref_cell = ref_grid_qua(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      if ( nodes[node_per] == faceid )
	{
	  fprintf(file," %d",nodes[4]);
	  fprintf(file,"\n");
	}

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    {
      node_per = ref_cell_node_per(ref_cell);
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  for ( node = 0; node < node_per; node++ )
	    fprintf(file," %d",o2n[nodes[node]]+1);
	  fprintf(file,"\n");
	}
    }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_export_bin_ugrid( REF_GRID ref_grid, char *filename,
					REF_BOOL swap)
{
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT nnode,ntri,nqua,ntet,npyr,npri,nhex;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell;
  REF_DBL swapped_dbl;
  REF_INT group;
  REF_INT faceid, min_faceid, max_faceid;

  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  nnode = ref_node_n(ref_node);

  ntri = ref_cell_n(ref_grid_tri(ref_grid));
  nqua = ref_cell_n(ref_grid_qua(ref_grid));

  ntet = ref_cell_n(ref_grid_tet(ref_grid));
  npyr = ref_cell_n(ref_grid_pyr(ref_grid));
  npri = ref_cell_n(ref_grid_pri(ref_grid));
  nhex = ref_cell_n(ref_grid_hex(ref_grid));

  if (swap) SWAP_INT(nnode);
  if (swap) SWAP_INT(ntri);
  if (swap) SWAP_INT(nqua);
  if (swap) SWAP_INT(ntet);
  if (swap) SWAP_INT(npyr);
  if (swap) SWAP_INT(npri);
  if (swap) SWAP_INT(nhex);

  REIS(1, fwrite(&nnode,sizeof(REF_INT),1,file),"nnode");

  REIS(1, fwrite(&ntri,sizeof(REF_INT),1,file),"ntri");
  REIS(1, fwrite(&nqua,sizeof(REF_INT),1,file),"nqua");

  REIS(1, fwrite(&ntet,sizeof(REF_INT),1,file),"ntet");
  REIS(1, fwrite(&npyr,sizeof(REF_INT),1,file),"npyr");
  REIS(1, fwrite(&npri,sizeof(REF_INT),1,file),"npri");
  REIS(1, fwrite(&nhex,sizeof(REF_INT),1,file),"nhex");

  RSS( ref_node_compact( ref_node, &o2n, &n2o), "compact" );

  for ( node = 0; node < ref_node_n(ref_node); node++ )
    {
      swapped_dbl = ref_node_xyz(ref_node,0,n2o[node]);
      if (swap) SWAP_DBL(swapped_dbl);
      REIS(1, fwrite(&swapped_dbl,sizeof(REF_DBL),1,file),"x");
      swapped_dbl = ref_node_xyz(ref_node,1,n2o[node]);
      if (swap) SWAP_DBL(swapped_dbl);
      REIS(1, fwrite(&swapped_dbl,sizeof(REF_DBL),1,file),"y");
      swapped_dbl = ref_node_xyz(ref_node,2,n2o[node]);
      if (swap) SWAP_DBL(swapped_dbl);
      REIS(1, fwrite(&swapped_dbl,sizeof(REF_DBL),1,file),"z");
    }

  RSS( ref_export_faceid_range( ref_grid, &min_faceid, &max_faceid), "range");

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      if ( nodes[node_per] == faceid )
	for ( node = 0; node < node_per; node++ )
	  {
	    nodes[node] = o2n[nodes[node]]+1;
	    if (swap) SWAP_INT(nodes[node]);
	    REIS(1, fwrite(&(nodes[node]),sizeof(REF_INT),1,file),"tri");
	  }

  ref_cell = ref_grid_qua(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      if ( nodes[node_per] == faceid )
	for ( node = 0; node < node_per; node++ )
	  {
	    nodes[node] = o2n[nodes[node]]+1;
	    if (swap) SWAP_INT(nodes[node]);
	    REIS(1, fwrite(&(nodes[node]),sizeof(REF_INT),1,file),"qua");
	  }

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      if ( nodes[node_per] == faceid )
	{
	  if (swap) SWAP_INT(nodes[3]);
	  REIS(1, fwrite(&(nodes[3]),sizeof(REF_INT),1,file),"tri id");
	}

  ref_cell = ref_grid_qua(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      if ( nodes[node_per] == faceid )
	{
	  if (swap) SWAP_INT(nodes[4]);
	  REIS(1, fwrite(&(nodes[4]),sizeof(REF_INT),1,file),"qua id");
	}

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    {
      node_per = ref_cell_node_per(ref_cell);
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	for ( node = 0; node < node_per; node++ )
	  {
	    nodes[node] = o2n[nodes[node]]+1;
	    if (swap) SWAP_INT(nodes[node]);
	    REIS(1, fwrite(&(nodes[node]),sizeof(REF_INT),1,file),"cell");
	  }
    }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_lb8_ugrid( REF_GRID ref_grid, char *filename  )
{
  RSS( ref_export_bin_ugrid( ref_grid, filename, REF_FALSE ), 
       "bin not swapped" );
  return REF_SUCCESS;
}

REF_STATUS ref_export_b8_ugrid( REF_GRID ref_grid, char *filename  )
{
  RSS( ref_export_bin_ugrid( ref_grid, filename, REF_TRUE ), 
       "bin swap" );
  return REF_SUCCESS;
}

REF_STATUS ref_export_faceid_range( REF_GRID ref_grid, 
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

  return REF_SUCCESS;
}

REF_STATUS ref_export_cogsg( REF_GRID ref_grid, char *filename_cogsg )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node, cell, ixyz;
  REF_INT nnode, nbn, ntri;
  REF_INT *o2n, *n2o;
  REF_CELL ref_cell;
  REF_DICT ref_dict;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  FILE *file;
  REF_INT boundary_index, boundary_tag;

  int fortran_record_size;
  REF_INT i;
  REF_DBL t;
  size_t end_of_string;
  char filename_bc[1024];
  char *period_bc;
  end_of_string = strlen(filename_cogsg);

  RAS( end_of_string > 6, "filename too short" );
  RAS( end_of_string < 1023, "filename too long" );
  REIS(0, strcmp(&filename_cogsg[end_of_string-6],".cogsg"), 
       "filename must end in .cogsg" );

  sprintf(filename_bc, "%s", filename_cogsg);
  period_bc = strrchr(filename_bc, '.');
  sprintf(period_bc, ".bc" );

  REIS( 0, ref_cell_n( ref_grid_qua(ref_grid) ), "no quad support");
  REIS( 0, ref_cell_n( ref_grid_pyr(ref_grid) ), "no pyramid support");
  REIS( 0, ref_cell_n( ref_grid_pri(ref_grid) ), "no prism support");
  REIS( 0, ref_cell_n( ref_grid_hex(ref_grid) ), "no hex support");

  ref_malloc_init( o2n, ref_node_max(ref_node), REF_INT, REF_EMPTY );
  ref_malloc( n2o, ref_node_n(ref_node), REF_INT );

  nnode = 0;

  ref_cell = ref_grid_tri(ref_grid);

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    for ( node = 0 ; node < ref_cell_node_per(ref_cell) ; node++ )
      if ( REF_EMPTY == o2n[nodes[node]] )
	{
	  o2n[nodes[node]] = nnode;
	  nnode++;
	}
  nbn = nnode;
  
  each_ref_node_valid_node( ref_node, node )
    if ( REF_EMPTY == o2n[node] )
      {
	o2n[node] = nnode;
	nnode++;
      }

  REIS( nnode, ref_node_n(ref_node), "nnode miscount" );

  each_ref_node_valid_node( ref_node, node )
    n2o[o2n[node]] = node;

  RSS( ref_dict_create( &ref_dict ), "create dict" ); 

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    RSS( ref_dict_store( ref_dict, nodes[3], REF_EMPTY ), "mark tri" );

  /*
  2)  project.bc: patch/surface-triangle file (ASCII)
       integer fapch(mbf),fnode(mbf,3)
       character*1 text(80)
       open(unit=12,file='project.bc',form='formatted')

       read(12,*)nbf,nbc,npatch,igrid
       read(12,900)text
       do if=1,nbf
        read(12,*)jf,fapch(if),fnode(if,in),in=1,3)
       enddo
  900  format(80a1)

        nbf = number of boundary triangular faces
        nbc = number of surface grid nodes along patch boundaries (curves)
     npatch = number of surface patches
      igrid = 1 for inviscid grids; 2 for viscous grids
       text = text line
         jf = triangle index
  fapch(if) = surface patch index containing surface triangle "if"
fnode(if,in) = node "in" of triangle "if"

Note: triangle connectivities are according to the right-hand rule with
      the outward normals pointing into the computational domain.
   */

  file = fopen(filename_bc,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename_bc);
  RNS(file, "unable to open file" );

  fprintf(file,"%d %d %d %d\n",
	  ref_cell_n(ref_cell),
	  0,
	  ref_dict_n(ref_dict),
	  1);
  fprintf(file,"exported by ref_export_cogsg x         x         x         x         x        80\n");
 
  ntri = 0;

  each_ref_dict_key( ref_dict, boundary_index, boundary_tag )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( boundary_tag == nodes[3] )
      {
	fprintf(file," %d %d",ntri+1,boundary_tag);
	ntri++;
	for ( node = 0; node < 3; node++ )
	  {
	    fprintf(file," %d",o2n[nodes[node]] + 1);
	  }
	fprintf(file,"\n");
      }

  fclose(file);

  RSS( ref_dict_free( ref_dict), "free");

  /*
  3)  project.cogsg: x,y,z coordinates of the grid nodes and
                    tetrahedral node connectivity
      real(8)  crd(mp,3),t
      integer int(mc,4)
       open(9,file='project.cogrd',form='unformatted',iostat=ios,
     &      err=555,status='old')

      read(9)inew,nc,npo,nbn,npv,nev,t,
     &       ((int(ie,in),ie=1,nc),in=1,4)
      read(9)((crd(ip,id),ip=1,npo),id=1,3)
      read(9)0
      read(9)0

  where

         inew = a dummy variable (integer) should be -1
           nc = number of tetrahedral cells
          npo = total number of grid nodes (including nbn)
          nbn = number of grid nodes on the boundaries (including nbc)
          npv = number of grid points in the viscous layers
                (=0 for Euler grids)
          ncv = number of cells in the viscous layers
                (=0 for Euler grids)
            t = a dummy variable (real - double)
   int(ie,in) = tetradhedral cell connectivity
                (node "in" of cell "ie")
   crd(ip,id) = x, y, and z coordinates of node "ip"

  Note 1: the first "nbn" coordinates listed in this file are those of
         the boundary nodes.

  Note 2:  tetrahedral cell connectivities are given according to the
           right-hand rule (3 nodes of the base in the counter-clockwise
           direction followed by the 4th node.)
   */

  file = fopen(filename_cogsg,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename_cogsg);
  RNS(file, "unable to open file" );

  ref_cell = ref_grid_tet(ref_grid);

  fortran_record_size = 4*6 + 8*1 + 4*4*ref_cell_n(ref_cell);

  i = fortran_record_size; 
  SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"record");
  
  i = -1; SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"inew");
  i = ref_cell_n(ref_cell); SWAP_INT(i); 
  REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"nc");
  i = nnode; SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"npo");
  i = nbn; SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"nbn");
  i = 0; SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"npv");
  i = 0; SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"nev");
  t = 0.0; SWAP_DBL(t); REIS(1, fwrite(&t,sizeof(REF_DBL),1,file),"t");

  for ( node = 0 ; node < ref_cell_node_per(ref_cell) ; node++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      {
	i = o2n[nodes[node]] + 1;
	SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"int");
      }

  i = fortran_record_size; 
  SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"record");

  fortran_record_size = 8*3*nnode;

  i = fortran_record_size; 
  SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"record");
  for ( ixyz = 0 ; ixyz < 3 ; ixyz++ )
    for ( node = 0 ; node < nnode ; node++ )
      {
	t = ref_node_xyz(ref_node,ixyz,n2o[node]);
	SWAP_DBL(t); REIS(1, fwrite(&t,sizeof(REF_DBL),1,file),"crd");
      }
  i = fortran_record_size; 
  SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"record");

  i = 4; SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"record");
  i = 0; SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"nelec");
  i = 4; SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"record");

  i = 4; SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"record");
  i = 0; SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"npoic");
  i = 4; SWAP_INT(i); REIS(1, fwrite(&i,sizeof(REF_INT),1,file),"record");

  fclose(file);

  ref_free( o2n );
  ref_free( n2o );

  return REF_SUCCESS;
}


REF_STATUS ref_export_c( REF_GRID ref_grid, char *filename  )
{
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node, ixyz;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;

  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  RSS( ref_node_compact( ref_node, &o2n, &n2o), "compact" );

  for ( node = 0; node < ref_node_n(ref_node); node++ )
    {
      fprintf(file, "  RSS(ref_node_add(ref_node,%d,&node),\"node\");\n", node);
      for ( ixyz = 0; ixyz < 3; ixyz++ ) 
	fprintf(file, "  ref_node_xyz(ref_node,%d,node) = %.15e;\n",
		ixyz, ref_node_xyz(ref_node,ixyz,n2o[node]));
    }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      for ( node = 0; node < ref_cell_size_per(ref_cell); node++ )
	fprintf(file,"  nodes[%d] = %d;\n",node,o2n[nodes[node]]);
      fprintf(file,"  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),\"qua\");\n");
    }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_pdf( REF_GRID ref_grid, char *pdf_filename )
{
  char temp_filename[] = "ref_export_temp_for_pdf.eps";
  char command[1024];
  RSS( ref_export_eps( ref_grid, temp_filename ), "temp eps");
  sprintf( command, "epstopdf %s -o=%s", temp_filename, pdf_filename );
  REIS(0, system( command ), "epstopdf failed");
  REIS(0, remove( temp_filename ), "temp clean up");

  return REF_SUCCESS;
}

REF_STATUS ref_export_eps( REF_GRID ref_grid, char *filename )
{
  FILE *f;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT node, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  f = fopen("ref_export.gnuplot","w");
  if (NULL == (void *)f) printf("unable to open ref_export.gnuplot\n");
  RNS(f, "unable to open file" );

  fprintf(f,"reset\n");
  fprintf(f,"set term postscript eps\n");
  fprintf(f,"set output '%s'\n",filename);
  fprintf(f,"set size ratio -1\n");
  fprintf(f,"set xlabel 'X'\n");
  fprintf(f,"set ylabel 'Z'\n");
  fprintf(f,"plot '-' title '' with lines lw 0.5\n");

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      for ( node = 0; node < ref_cell_node_per(ref_cell); node++ )
	{
	  fprintf(f,"%25.15f %25.15f\n", 
		  ref_node_xyz(ref_node,0,nodes[node]),
		  ref_node_xyz(ref_node,2,nodes[node]));
	}
      node = 0;
      fprintf(f,"%25.15f %25.15f\n", 
	      ref_node_xyz(ref_node,0,nodes[node]),
	      ref_node_xyz(ref_node,2,nodes[node]));
      fprintf(f,"\n\n");

    }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      for ( node = 0; node < ref_cell_node_per(ref_cell); node++ )
	{
	  fprintf(f,"%25.15f %25.15f\n", 
		  ref_node_xyz(ref_node,0,nodes[node]),
		  ref_node_xyz(ref_node,2,nodes[node]));
	}
      node = 0;
      fprintf(f,"%25.15f %25.15f\n", 
	      ref_node_xyz(ref_node,0,nodes[node]),
	      ref_node_xyz(ref_node,2,nodes[node]));
      fprintf(f,"\n\n");

    }

  fprintf(f,"e\n");
  fclose(f);

  REIS(0, system("gnuplot ref_export.gnuplot"),"gnuplot failed");
  REIS(0, remove( "ref_export.gnuplot" ), "temp clean up");

  return REF_SUCCESS;
}

REF_STATUS ref_export_html( REF_GRID ref_grid, char *filename )
{
  FILE *f;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT node, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT *o2n, *n2o;

  f = fopen(filename,"w");
  if (NULL == (void *)f) printf("unable to open %s\n",filename);
  RNS(f, "unable to open file" );

  RSS( ref_node_compact( ref_node, &o2n, &n2o), "compact" );

  fprintf(f,"<html>\n");
  fprintf(f,"  <head>\n");
  fprintf(f,"    <title>refine export</title>\n");
  fprintf(f,"    <link rel='stylesheet' type='text/css'\n");
  fprintf(f,"          href='http://www.x3dom.org/download/x3dom.css'>\n");
  fprintf(f,"    </link>\n");
  fprintf(f,"    <script type='text/javascript'\n");
  fprintf(f,"          src='http://www.x3dom.org/download/x3dom.js'>\n");
  fprintf(f,"    </script>\n");
  fprintf(f,"    <style>\n");
  fprintf(f,"      x3d {width:100%%;height:100%%;border:none}\n");
  fprintf(f,"      body {margin:0;width:100%%;height:100%%;}\n");
  fprintf(f,"    </style>\n");
  fprintf(f,"  </head>\n");
  fprintf(f,"  <body id='body'>\n");
  fprintf(f,"    <a href=\"http://x3dom.org/docs/dev/navigation.html\">\n");
  fprintf(f,"camera control help\n");
  fprintf(f,"    </a>\n");
  fprintf(f,"    <x3d id='x3d'><scene><shape>\n");

  fprintf(f,"      <IndexedLineSet coordIndex='\n");
  if ( REF_TRUE )
    {
      ref_cell = ref_grid_tri(ref_grid);
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  for ( node = 0; node < ref_cell_node_per(ref_cell); node++ )
	    fprintf(f," %d",o2n[nodes[node]]);
	  fprintf(f," %d %d\n",o2n[nodes[0]],-1);
	}
    }
     
  if ( REF_TRUE )
    {
      ref_cell = ref_grid_qua(ref_grid);
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  for ( node = 0; node < ref_cell_node_per(ref_cell); node++ )
	    fprintf(f," %d",o2n[nodes[node]]);
	  fprintf(f," %d %d\n",o2n[nodes[0]],-1);
	}
    }

  if ( REF_TRUE )
    {
      ref_cell = ref_grid_tet(ref_grid);
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  fprintf(f," %d %d %d %d %d %d %d %d\n",
		  o2n[nodes[0]],
		  o2n[nodes[1]],
		  o2n[nodes[2]],
		  o2n[nodes[3]],
		  o2n[nodes[0]],
		  o2n[nodes[3]],
		  o2n[nodes[1]],
		  -1);
	}
    }

  fprintf(f,"      ' >\n");

  fprintf(f,"      <Coordinate point='\n");
  for ( node = 0; node < ref_node_n(ref_node); node++ )
    {
      fprintf(f," %.15e",ref_node_xyz(ref_node,0,n2o[node]));
      fprintf(f," %.15e",ref_node_xyz(ref_node,1,n2o[node]));
      fprintf(f," %.15e\n",ref_node_xyz(ref_node,2,n2o[node]));
    }
  fprintf(f,"      ' />\n");

  fprintf(f,"      </IndexedLineSet>\n");

  fprintf(f,"    </shape></scene></x3d>\n");
  fprintf(f,"  </body>\n");
  fprintf(f,"</html>\n");

  ref_free(n2o);
  ref_free(o2n);

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_export_meshb( REF_GRID ref_grid, char *filename )
{
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT *o2n, *n2o;
  REF_INT code, version, keyword_code, dim;
  int next_position;
  REF_INT node;
  REF_INT min_faceid, max_faceid, node_per, faceid, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT id;

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  RSS( ref_node_compact( ref_node, &o2n, &n2o), "compact" );

  code = 1;
  REIS(1, fwrite(&code,sizeof(int),1,file),"code");
  version = 2;
  REIS(1, fwrite(&version,sizeof(int),1,file),"version");
  next_position = 4+4+4+ftell(file);
  keyword_code = 3;
  REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"dim code");
  REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");
  dim = 3;
  REIS(1, fwrite(&dim,sizeof(int),1,file),"dim");

  next_position = 4+4+4+ref_node_n(ref_node)*(3*8+4)+ftell(file);
  keyword_code = 4;
  REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"vertex version code");
  REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");
  REIS(1, fwrite(&(ref_node_n(ref_node)),sizeof(int),1,file),"nnode");

  for ( node = 0; node < ref_node_n(ref_node); node++ )
    {
      REIS(1, fwrite(&(ref_node_xyz(ref_node,0,n2o[node])),
		     sizeof(double),1,file),"x");
      REIS(1, fwrite(&(ref_node_xyz(ref_node,1,n2o[node])),
		     sizeof(double),1,file),"y");
      REIS(1, fwrite(&(ref_node_xyz(ref_node,2,n2o[node])),
		     sizeof(double),1,file),"z");
      id = node+1;
      REIS(1, fwrite(&(id),sizeof(int),1,file),"id");
    }

  ref_cell = ref_grid_tri(ref_grid);

  next_position = 4+4+4+ref_cell_n(ref_cell)*(4*4)+ftell(file);
  keyword_code = 6;
  REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"vertex version code");
  REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");
  REIS(1, fwrite(&(ref_cell_n(ref_cell)),sizeof(int),1,file),"nnode");

  RSS( ref_export_faceid_range( ref_grid, &min_faceid, &max_faceid), "range");

  node_per = ref_cell_node_per(ref_cell);
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      if ( nodes[node_per] == faceid )
	{
	  for ( node = 0; node < node_per; node++ )
	    {
	      nodes[node] = o2n[nodes[node]]+1;
	      REIS(1, fwrite(&(nodes[node]),sizeof(REF_INT),1,file),"tri");
	    }
	  REIS(1, fwrite(&(nodes[3]),sizeof(REF_INT),1,file),"tri id");
	}

  ref_cell = ref_grid_tet(ref_grid);

  next_position = 4+4+4+ref_cell_n(ref_cell)*(4*5)+ftell(file);
  keyword_code = 8;
  REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"vertex version code");
  REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");
  REIS(1, fwrite(&(ref_cell_n(ref_cell)),sizeof(int),1,file),"nnode");
  node_per = ref_cell_node_per(ref_cell);
  id = 0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      id++;
      for ( node = 0; node < node_per; node++ )
	{
	  nodes[node] = o2n[nodes[node]]+1;
	  REIS(1, fwrite(&(nodes[node]),sizeof(int),1,file),"cell");
	}
      REIS(1, fwrite(&(id),sizeof(int),1,file),"tri id");
    }

  /* End */
  keyword_code = 54;
  REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"vertex version code");
  next_position = 0;
  REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}
REF_STATUS ref_export_twod_meshb( REF_GRID ref_grid, char *filename )
{
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT nnode, nedge, ntri;
  REF_BOOL twod_node, twod_edge;
  REF_INT *o2n, *n2o;
  REF_INT code, version, keyword_code, dim;
  int next_position;
  REF_INT node;
  REF_INT min_faceid, max_faceid, node_per, faceid, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT id;
  REF_INT node0, node1, node0_1, node1_1;
  
  RAS( ref_grid_twod(ref_grid), "expected twod convention grid" );
  
  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  ref_malloc_init( o2n, ref_node_max(ref_node), REF_INT, REF_EMPTY );
  ref_malloc_init( n2o, ref_node_max(ref_node), REF_INT, REF_EMPTY );

  nnode = 0;
  each_ref_node_valid_node( ref_node, node )
    {
      RSS( ref_node_node_twod( ref_node, node, &twod_node ), "twod node" );
      if ( twod_node ) 
	{
	  o2n[node] = nnode;
	  n2o[nnode] = node;
	  nnode++;
	}
    }

  code = 1;
  REIS(1, fwrite(&code,sizeof(int),1,file),"code");
  version = 2;
  REIS(1, fwrite(&version,sizeof(int),1,file),"version");
  next_position = 4+4+4+ftell(file);
  keyword_code = 3;
  REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"dim code");
  REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");
  dim = 2;
  REIS(1, fwrite(&dim,sizeof(int),1,file),"dim");

  next_position = 4+4+4+nnode*(2*8+4)+ftell(file);
  keyword_code = 4;
  REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"vertex version code");
  REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");
  REIS(1, fwrite(&(nnode),sizeof(int),1,file),"nnode");

  for ( node = 0; node < nnode; node++ )
    {
      REIS(1, fwrite(&(ref_node_xyz(ref_node,0,n2o[node])),
		     sizeof(double),1,file),"x");
      REIS(1, fwrite(&(ref_node_xyz(ref_node,2,n2o[node])),
		     sizeof(double),1,file),"z");
      id = node+1;
      REIS(1, fwrite(&(id),sizeof(int),1,file),"id");
    }

  ref_cell = ref_grid_qua(ref_grid);

  next_position = 4+4+4+ref_cell_n(ref_cell)*(3*4)+ftell(file);
  keyword_code = 5;
  REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"vertex version code");
  REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");
  REIS(1, fwrite(&(ref_cell_n(ref_cell)),sizeof(int),1,file),"nnode");

  RSS( ref_export_faceid_range( ref_grid, &min_faceid, &max_faceid), "range");

  nedge=0;
  node_per = ref_cell_node_per(ref_cell);
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      if ( nodes[node_per] == faceid )
	{
	  for ( node0 = 0; node0 < node_per; node0++ )
	    {

	      node1 = node0+1;
	      if ( node1>3 ) node1 = 0;
	      RSS( ref_node_edge_twod( ref_node, 
				       nodes[node0], 
				       nodes[node1],
				       &twod_edge ), "twod edge" );
	      if ( twod_edge )
		{
		  nedge++;
		  node0_1 = o2n[nodes[node0]]+1;
		  node1_1 = o2n[nodes[node1]]+1;
		  REIS(1, fwrite(&(node0_1),sizeof(REF_INT),1,file),"edge n0");
		  REIS(1, fwrite(&(node1_1),sizeof(REF_INT),1,file),"edge n1");
		}
	    }
	  REIS(1, fwrite(&(nodes[4]),sizeof(REF_INT),1,file),"edge id");
	}
  REIS( nedge, ref_cell_n(ref_cell), "edge/quad miscount" );

  ref_cell = ref_grid_tri(ref_grid);
  ntri = ref_cell_n(ref_cell)/2;
  
  next_position = 4+4+4+ntri*(4*4)+ftell(file);
  keyword_code = 6;
  REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"vertex version code");
  REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");
  REIS(1, fwrite(&(ntri),sizeof(int),1,file),"nnode");

  RSS( ref_export_faceid_range( ref_grid, &min_faceid, &max_faceid), "range");

  ntri = 0;
  node_per = ref_cell_node_per(ref_cell);
  for ( faceid = min_faceid ; faceid <= max_faceid ; faceid++ )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      if ( nodes[node_per] == faceid )
	{
	  RSS( ref_node_node_twod( ref_node, nodes[0], &twod_node ), "2d nod" );
	  if ( twod_node )
	    {
	      ntri++;
	      for ( node = 0; node < node_per; node++ )
		{
		  nodes[node] = o2n[nodes[node]]+1;
		  REIS(1, fwrite(&(nodes[node]),sizeof(REF_INT),1,file),"tri");
		}
	      REIS(1, fwrite(&(nodes[3]),sizeof(REF_INT),1,file),"tri id");
	    }
	}
  REIS( ntri, ref_cell_n(ref_cell)/2, "triangle miscount" );
  
  /* End */
  keyword_code = 54;
  REIS(1, fwrite(&keyword_code,sizeof(int),1,file),"vertex version code");
  next_position = 0;
  REIS(1, fwrite(&next_position,sizeof(next_position),1,file),"next pos");

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_twod_msh( REF_GRID ref_grid, char *filename )
{
  FILE *f;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_BOOL twod_node;
  REF_INT nnode;
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node0, node1;
  REF_INT nedge;
  REF_BOOL twod_edge;
  REF_INT ntri;

  RAS( ref_grid_twod(ref_grid), "expected twod convention grid" );

  f = fopen(filename,"w");
  if (NULL == (void *)f) printf("unable to open %s\n",filename);
  RNS(f, "unable to open file" );

  fprintf(f, "MeshVersionFormatted 0\n\n");
  fprintf(f, "Dimension 2\n\n");

  ref_malloc_init( o2n, ref_node_max(ref_node), REF_INT, REF_EMPTY );
  ref_malloc_init( n2o, ref_node_max(ref_node), REF_INT, REF_EMPTY );
 
  nnode = 0;
  each_ref_node_valid_node( ref_node, node )
    {
      RSS( ref_node_node_twod( ref_node, node, &twod_node ), "twod node" );
      if ( twod_node ) 
	{
	  o2n[node] = nnode;
	  n2o[nnode] = node;
	  nnode++;
	}
    }

  fprintf(f, "\nVertices\n%d\n", nnode );
  for ( node = 0; node < nnode; node++ )
    {
      fprintf(f, "%.16E %.16E %d\n", 
	      ref_node_xyz(ref_node,0,n2o[node]), 
	      ref_node_xyz(ref_node,2,n2o[node]), 
	      1);
    }

  ref_cell = ref_grid_qua(ref_grid);
  fprintf(f, "\nEdges\n%d\n", ref_cell_n(ref_cell) );
  nedge=0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      for ( node0 = 0; node0 < 4; node0++ )
	{
	  node1 = node0+1;
	  if ( node1>3 ) node1 = 0;
	  RSS( ref_node_edge_twod( ref_node, 
				   nodes[node0], 
				   nodes[node1],
				   &twod_edge ), "twod edge" );
	  if ( twod_edge )
	    {
	      nedge++;
	      fprintf(f, "%d %d %d\n", 
		      o2n[nodes[node0]]+1, 
		      o2n[nodes[node1]]+1, 
		      nodes[4]);
	    }
	}
    }
  REIS( nedge, ref_cell_n(ref_cell), "edge/quad miscount" );

  ref_cell = ref_grid_tri(ref_grid);
  fprintf(f, "\nTriangles\n%d\n", ref_cell_n(ref_cell)/2 );
  ntri = 0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      RSS( ref_node_node_twod( ref_node, nodes[0], &twod_node ), "twod node" );
      if ( twod_node )
	{
	  ntri++;
	  fprintf(f, "%d %d %d %d\n", 
		  o2n[nodes[0]]+1, 
		  o2n[nodes[2]]+1, 
		  o2n[nodes[1]]+1, 
		  nodes[3]);
	}
    }
  REIS( ntri, ref_cell_n(ref_cell)/2, "triangle miscount" );

  ref_free(n2o);
  ref_free(o2n);

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_export_metric2d( REF_GRID ref_grid, char *filename )
{
  FILE *f;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_BOOL twod_node;
  REF_INT nnode;

  f = fopen(filename,"w");
  if (NULL == (void *)f) printf("unable to open %s\n",filename);
  RNS(f, "unable to open file" );

  ref_malloc_init( o2n, ref_node_max(ref_node), REF_INT, REF_EMPTY );
  ref_malloc_init( n2o, ref_node_max(ref_node), REF_INT, REF_EMPTY );
 
  nnode = 0;
  each_ref_node_valid_node( ref_node, node )
    {
      RSS( ref_node_node_twod( ref_node, node, &twod_node ), "twod node" );
      if ( twod_node ) 
	{
	  o2n[node] = nnode;
	  n2o[nnode] = node;
	  nnode++;
	}
    }

  fprintf(f, "%d %d\n", nnode, 3 );
  for ( node = 0; node < nnode; node++ )
    {
      fprintf(f, "%.16E %.16E  %.16E \n", 
	      ref_node_metric(ref_node,0,n2o[node]), 
	      ref_node_metric(ref_node,2,n2o[node]), 
	      ref_node_metric(ref_node,5,n2o[node]) );
    }

  ref_free(n2o);
  ref_free(o2n);

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_export_twod_sol( REF_GRID ref_grid, char *filename )
{
  FILE *f;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_BOOL twod_node;
  REF_INT nnode;

  f = fopen(filename,"w");
  if (NULL == (void *)f) printf("unable to open %s\n",filename);
  RNS(f, "unable to open file" );

  ref_malloc_init( o2n, ref_node_max(ref_node), REF_INT, REF_EMPTY );
  ref_malloc_init( n2o, ref_node_max(ref_node), REF_INT, REF_EMPTY );
 
  nnode = 0;
  each_ref_node_valid_node( ref_node, node )
    {
      RSS( ref_node_node_twod( ref_node, node, &twod_node ), "twod node" );
      if ( twod_node ) 
	{
	  o2n[node] = nnode;
	  n2o[nnode] = node;
	  nnode++;
	}
    }

  fprintf(f,"MeshVersionFormatted 2\n\n");
  fprintf(f,"Dimension 2\n\n");
  fprintf(f,"SolAtVertices\n%d\n1 3\n",nnode);
  
  for ( node = 0; node < nnode; node++ )
    {
      fprintf(f, "%.16E %.16E  %.16E \n", 
	      ref_node_metric(ref_node,0,n2o[node]), 
	      ref_node_metric(ref_node,2,n2o[node]), 
	      ref_node_metric(ref_node,5,n2o[node]) );
    }

  ref_free(n2o);
  ref_free(o2n);

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_export_plt( REF_GRID ref_grid, char *filename  )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  FILE *file;
  int one = 1;
  int filetype = 0;
  int ascii[8];
  int numvar = 3;
  float zonemarker = 299.0;
  int parentzone = -1;
  int strandid = -1;
  double solutiontime = 0.0;
  int notused = -1;
  int zonetype = 4; /*4=FETETRAHEDRON*/
  int datapacking = 1; /*1=Point*/
  int varloc = 0; /*0 = Don't specify, all data is located at nodes*/
  int faceneighbors = 0;
  int numpts = ref_node_n(ref_node);
  int numelements = ref_cell_n(ref_cell);
  int celldim = 0;
  int aux = 0;
  float eohmarker = 357.0;
  int dataformat = 1;
  int passive = 0;
  int varsharing = 0;
  int connsharing = -1;
  float data;

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  printf("%p\n",(void *)ref_grid);

  REIS(8, fwrite(&"#!TDV112",sizeof(char),8,file),"header");
  REIS(1, fwrite(&one,sizeof(int),1,file),"magic");
  REIS(1, fwrite(&filetype,sizeof(int),1,file),"filetype");

  ascii[0] = (int)'f';
  ascii[1] = (int)'t';
  ascii[2] = 0;
  REIS(3, fwrite(&ascii,sizeof(int),3,file),"title");

  REIS(1, fwrite(&numvar,sizeof(int),1,file),"numvar");
  ascii[0] = (int)'x';
  ascii[1] = 0;
  REIS(2, fwrite(&ascii,sizeof(int),2,file),"var");
  ascii[0] = (int)'y';
  ascii[1] = 0;
  REIS(2, fwrite(&ascii,sizeof(int),2,file),"var");
  ascii[0] = (int)'z';
  ascii[1] = 0;
  REIS(2, fwrite(&ascii,sizeof(int),2,file),"var");

  REIS(1, fwrite(&zonemarker,sizeof(float),1,file),"zonemarker");

  ascii[0] = (int)'e';
  ascii[1] = (int)'4';
  ascii[2] = 0;
  REIS(3, fwrite(&ascii,sizeof(int),3,file),"title");

  REIS(1, fwrite(&parentzone,sizeof(int),1,file),"int");
  REIS(1, fwrite(&strandid,sizeof(int),1,file),"int");
  REIS(1, fwrite(&solutiontime,sizeof(double),1,file),"double");
  REIS(1, fwrite(&notused,sizeof(int),1,file),"int");
  REIS(1, fwrite(&zonetype,sizeof(int),1,file),"int");
  REIS(1, fwrite(&datapacking,sizeof(int),1,file),"int");
  REIS(1, fwrite(&varloc,sizeof(int),1,file),"int");
  REIS(1, fwrite(&faceneighbors,sizeof(int),1,file),"int");
  REIS(1, fwrite(&numpts,sizeof(int),1,file),"int");
  REIS(1, fwrite(&numelements,sizeof(int),1,file),"int");
  REIS(1, fwrite(&celldim,sizeof(int),1,file),"int");
  REIS(1, fwrite(&celldim,sizeof(int),1,file),"int");
  REIS(1, fwrite(&celldim,sizeof(int),1,file),"int");
  REIS(1, fwrite(&aux,sizeof(int),1,file),"int");

  REIS(1, fwrite(&eohmarker,sizeof(float),1,file),"eohmarker");
  REIS(1, fwrite(&zonemarker,sizeof(float),1,file),"zonemarker");

  REIS(1, fwrite(&dataformat,sizeof(int),1,file),"int");
  REIS(1, fwrite(&dataformat,sizeof(int),1,file),"int");
  REIS(1, fwrite(&dataformat,sizeof(int),1,file),"int");

  REIS(1, fwrite(&passive,sizeof(int),1,file),"int");
  REIS(1, fwrite(&varsharing,sizeof(int),1,file),"int");
  REIS(1, fwrite(&connsharing,sizeof(int),1,file),"int");

  fclose(file);
  return REF_SUCCESS;
}
