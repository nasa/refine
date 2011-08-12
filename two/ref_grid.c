
#include <stdlib.h>
#include <stdio.h>

#include "ref_grid.h"

REF_STATUS ref_grid_create( REF_GRID *ref_grid_ptr )
{
  (*ref_grid_ptr) = NULL;
  (*ref_grid_ptr) = (REF_GRID)malloc( sizeof(REF_GRID_STRUCT) );
  RNS(*ref_grid_ptr,"malloc ref_grid NULL");

  RSS( ref_node_create( &(*ref_grid_ptr)->nodes), "node create" );

  (*ref_grid_ptr)->cells[0] = NULL;
  (*ref_grid_ptr)->cells[1] = NULL;
  (*ref_grid_ptr)->cells[2] = NULL;
  (*ref_grid_ptr)->cells[3] = NULL;
  RSS( ref_cell_create( 4, &((*ref_grid_ptr)->cells[4]) ), "tet create" );
  RSS( ref_cell_create( 5, &((*ref_grid_ptr)->cells[5]) ), "pri create" );
  RSS( ref_cell_create( 6, &((*ref_grid_ptr)->cells[6]) ), "pyr create" );
  (*ref_grid_ptr)->cells[7] = NULL;
  RSS( ref_cell_create( 8, &((*ref_grid_ptr)->cells[8]) ), "hex create" );

  (*ref_grid_ptr)->faces[0] = NULL;
  (*ref_grid_ptr)->faces[1] = NULL;
  (*ref_grid_ptr)->faces[2] = NULL;
  RSS( ref_cell_create( 4, &((*ref_grid_ptr)->faces[3]) ), "tri create" );
  RSS( ref_cell_create( 5, &((*ref_grid_ptr)->faces[4]) ), "qua create" );

  return REF_SUCCESS;
}

REF_STATUS ref_grid_free( REF_GRID ref_grid )
{
  if ( NULL == (void *)ref_grid ) return REF_NULL;

  RSS( ref_node_free( ref_grid->nodes ), "node free");

  RFS( ref_cell_free( ref_grid->cells[0] ), "cell 0 free");
  RFS( ref_cell_free( ref_grid->cells[1] ), "cell 1 free");
  RFS( ref_cell_free( ref_grid->cells[2] ), "cell 2 free");
  RFS( ref_cell_free( ref_grid->cells[3] ), "cell 3 free");
  RSS( ref_cell_free( ref_grid->cells[4] ), "cell 4 free");
  RSS( ref_cell_free( ref_grid->cells[5] ), "cell 5 free");
  RSS( ref_cell_free( ref_grid->cells[6] ), "cell 6 free");
  RFS( ref_cell_free( ref_grid->cells[7] ), "cell 7 free");
  RSS( ref_cell_free( ref_grid->cells[8] ), "cell 8 free");

  RFS( ref_cell_free( ref_grid->faces[0] ), "face 0 free");
  RFS( ref_cell_free( ref_grid->faces[1] ), "face 1 free");
  RFS( ref_cell_free( ref_grid->faces[2] ), "face 2 free");
  RSS( ref_cell_free( ref_grid->faces[3] ), "face 3 free");
  RSS( ref_cell_free( ref_grid->faces[4] ), "face 4 free");

  ref_cond_free( ref_grid );
  return REF_SUCCESS;
}

REF_STATUS ref_grid_from_ugrid( char *filename, REF_GRID *ref_grid_ptr )
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
  ref_node = ref_grid_node(ref_grid);

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
