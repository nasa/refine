
#include <stdlib.h>
#include <stdio.h>

#include "ref_grid_import.h"

REF_STATUS ref_grid_import_ugrid( char *filename, REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_INT node, new_node;
  REF_DBL xyz[3];
  REF_INT tri, new_tri;
  REF_INT nodes[8];
  REF_INT qua, new_qua;
  REF_INT face_id;
  REF_INT cell, new_cell;

  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid->nodes;

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s",filename);
  RNS(file, "unable to open file" );

  RES( 1, fscanf( file, "%d", &nnode ), "nnode" );
  RES( 1, fscanf( file, "%d", &ntri ), "ntri" );
  RES( 1, fscanf( file, "%d", &nqua ), "nqua" );
  RES( 1, fscanf( file, "%d", &ntet ), "ntet" );
  RES( 1, fscanf( file, "%d", &npyr ), "npyr" );
  RES( 1, fscanf( file, "%d", &npri ), "npri" );
  RES( 1, fscanf( file, "%d", &nhex ), "nhex" );

  for( node=0; node<nnode ; node++ ) 
    {
      RSS( ref_node_add(ref_node, node, &new_node ), "new_node");
      RES( 1, fscanf( file, "%lf", &(xyz[0]) ), "x" );
      RES( 1, fscanf( file, "%lf", &(xyz[1]) ), "y" );
      RES( 1, fscanf( file, "%lf", &(xyz[2]) ), "z" );
      ref_node_xyz( ref_node, 0, new_node ) = xyz[0];
      ref_node_xyz( ref_node, 1, new_node ) = xyz[1];
      ref_node_xyz( ref_node, 2, new_node ) = xyz[2];
    }

  ref_cell = ref_grid->faces[3];
  nodes[3] = REF_EMPTY;
  for( tri = 0; tri < ntri ; tri++ )
    {
      for ( node = 0 ; node < 3 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "tri" );
      nodes[0]--; nodes[1]--; nodes[2]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_tri ), "new tri");
      RES( tri, new_tri, "tri index");
    }

  ref_cell = ref_grid->faces[4];
  nodes[4] = REF_EMPTY;
  for( qua = 0; qua < nqua ; qua++ )
    {
      for ( node = 0 ; node < 4 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "tri" );
      nodes[0]--; nodes[1]--; nodes[2]--; nodes[3]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_qua ), "new qua");
      RES( qua, new_qua, "qua index");
    }
  
  ref_cell = ref_grid->faces[3];
  for( tri = 0; tri < ntri ; tri++ )
    {
      RES( 1, fscanf( file, "%d", &face_id ), "tri id" );
      ref_cell_c2n(ref_cell,4,tri) = face_id;
    }

  ref_cell = ref_grid->faces[4];
  for( qua = 0; qua < nqua ; qua++ )
    {
      RES( 1, fscanf( file, "%d", &face_id ), "qua id" );
      ref_cell_c2n(ref_cell,5,qua) = face_id;
    }

  ref_cell = ref_grid->cells[4];
  for( cell = 0; cell < ntet ; cell++ )
    {
      for ( node = 0 ; node < 4 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "tet" );
      nodes[0]--; nodes[1]--; nodes[2]--; nodes[3]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "new tet");
      RES( cell, new_cell, "tet index");
    }

  ref_cell = ref_grid->cells[5];
  for( cell = 0; cell < npyr ; cell++ )
    {
      for ( node = 0 ; node < 5 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "pyr" );
      nodes[0]--; nodes[1]--; nodes[2]--; nodes[3]--; nodes[4]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "new pyr");
      RES( cell, new_cell, "pyr index");
    }

  ref_cell = ref_grid->cells[6];
  for( cell = 0; cell < npri ; cell++ )
    {
      for ( node = 0 ; node < 6 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "pri" );
      nodes[0]--; nodes[1]--; nodes[2]--; nodes[3]--; nodes[4]--; nodes[5]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "new pri");

      if ( cell != new_cell)
	{
	  printf("cell %d %d\n",cell,new_cell);
	}

      RES( cell, new_cell, "pri index");
    }

  ref_cell = ref_grid->cells[8];
  for( cell = 0; cell < nhex ; cell++ )
    {
      for ( node = 0 ; node < 8 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "hex" );
      nodes[0]--; nodes[1]--; nodes[2]--; nodes[3]--;
      nodes[4]--; nodes[5]--; nodes[6]--; nodes[7]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "new hex");
      RES( cell, new_cell, "hex index");
    }

  fclose(file);

  return REF_SUCCESS;
}
