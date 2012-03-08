
#include <stdlib.h>
#include <stdio.h>

#include <string.h>

#include "ref_import.h"

#include "ref_endian.h"

REF_STATUS ref_import_examine_header( char *filename )
{
  FILE *file;
  int i4, i4_swapped;
  long i8, i8_swapped;
  int i;

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  for(i=0;i<7;i++)
    {
      RES( 1, fread( &i4, sizeof(int), 1, file ), "int" );
      i4_swapped = i4; SWAP_INT(i4_swapped);
      printf(" %d: %d (%d swapped) ints\n",i,i4,i4_swapped);
    }

  rewind(file);

  for(i=0;i<3;i++)
    {
      RES( 1, fread( &i4, sizeof(int), 1, file ), "int" );
    }
  for(i=3;i<7;i++)
    {
      RES( 1, fread( &i8, sizeof(long), 1, file ), "long" );
      i8_swapped = i8; SWAP_INT(i8_swapped);
      printf(" %d: %ld (%ld swapped) long\n",i,i8,i8_swapped);
    }

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_import_by_extension( REF_GRID *ref_grid_ptr, char *filename )
{
  size_t end_of_string;

  end_of_string = strlen(filename);

  if( strcmp(&filename[end_of_string-9],".b8.ugrid") == 0 ) 
    {
      RSS( ref_import_b8_ugrid( ref_grid_ptr, filename ), "b8_ugrid failed");
    } 
  else 
    if( strcmp(&filename[end_of_string-6],".ugrid") == 0 ) 
      {
	RSS( ref_import_ugrid( ref_grid_ptr, filename ), "ugrid failed");
      } 
    else 
      if( strcmp(&filename[end_of_string-6],".fgrid") == 0 ) 
	{
	  RSS( ref_import_fgrid( ref_grid_ptr, filename ), "fgrid failed");
	} 
      else 
	if( strcmp(&filename[end_of_string-4],".msh") == 0 ) 
	  {
	    RSS( ref_import_msh( ref_grid_ptr, filename ), "msh failed");
	  } 
      else 
	{
	  printf("%s: %d: %s %s\n",__FILE__,__LINE__,
		 "input file name extension unknown", filename);
	  RSS( REF_FAILURE, "unknown file extension");
	}

  return REF_SUCCESS;
}

REF_STATUS ref_import_fgrid( REF_GRID *ref_grid_ptr, char *filename )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  FILE *file;
  REF_INT nnode, ntri, ntet;
  REF_INT ixyz, node, new_node;
  REF_DBL xyz;
  REF_INT tri, new_tri;
  REF_INT nodes[4];
  REF_INT face_id;
  REF_INT cell, new_cell;

  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  RES( 1, fscanf( file, "%d", &nnode ), "nnode" );
  RES( 1, fscanf( file, "%d", &ntri ), "ntri" );
  RES( 1, fscanf( file, "%d", &ntet ), "ntet" );

  for( node=0; node<nnode ; node++ )
    {
      RSS( ref_node_add(ref_node, node, &new_node ), "new_node");
      RES( node, new_node, "node index");
    }

  RSS( ref_node_initialize_n_global( ref_node, nnode ), "init glob");

  for( ixyz=0; ixyz<3 ; ixyz++ )
    for( node=0; node<nnode ; node++ ) 
      {
	RES( 1, fscanf( file, "%lf", &xyz ), "xyz" );
	ref_node_xyz( ref_node, ixyz, node ) = xyz;
      }

  ref_cell = ref_grid_tri(ref_grid);
  nodes[3] = REF_EMPTY;
  for( tri = 0; tri < ntri ; tri++ )
    {
      for ( node = 0 ; node < 3 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "tri" );
      nodes[0]--; nodes[1]--; nodes[2]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_tri ), "new tri");
      RES( tri, new_tri, "tri index");
    }

  ref_cell = ref_grid_tri(ref_grid);
  for( tri = 0; tri < ntri ; tri++ )
    {
      RES( 1, fscanf( file, "%d", &face_id ), "tri id" );
      ref_cell_c2n(ref_cell,3,tri) = face_id;
    }

  ref_cell = ref_grid_tet(ref_grid);
  for( cell = 0; cell < ntet ; cell++ )
    {
      for ( node = 0 ; node < 4 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "tet" );
      nodes[0]--; nodes[1]--; nodes[2]--; nodes[3]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "new tet");
      RES( cell, new_cell, "tet index");
    }

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_import_ugrid( REF_GRID *ref_grid_ptr, char *filename )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_INT node, new_node;
  REF_DBL xyz[3];
  REF_INT tri, new_tri;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT qua, new_qua;
  REF_INT face_id;
  REF_INT cell, new_cell;

  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
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
      RES( node, new_node, "node index");
      RES( 1, fscanf( file, "%lf", &(xyz[0]) ), "x" );
      RES( 1, fscanf( file, "%lf", &(xyz[1]) ), "y" );
      RES( 1, fscanf( file, "%lf", &(xyz[2]) ), "z" );
      ref_node_xyz( ref_node, 0, new_node ) = xyz[0];
      ref_node_xyz( ref_node, 1, new_node ) = xyz[1];
      ref_node_xyz( ref_node, 2, new_node ) = xyz[2];
    }

  RSS( ref_node_initialize_n_global( ref_node, nnode ), "init glob");

  ref_cell = ref_grid_tri(ref_grid);
  nodes[3] = REF_EMPTY;
  for( tri = 0; tri < ntri ; tri++ )
    {
      for ( node = 0 ; node < 3 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "tri" );
      nodes[0]--; nodes[1]--; nodes[2]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_tri ), "new tri");
      RES( tri, new_tri, "tri index");
    }

  ref_cell = ref_grid_qua(ref_grid);
  nodes[4] = REF_EMPTY;
  for( qua = 0; qua < nqua ; qua++ )
    {
      for ( node = 0 ; node < 4 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "qua" );
      nodes[0]--; nodes[1]--; nodes[2]--; nodes[3]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_qua ), "new qua");
      RES( qua, new_qua, "qua index");
    }
  
  ref_cell = ref_grid_tri(ref_grid);
  for( tri = 0; tri < ntri ; tri++ )
    {
      RES( 1, fscanf( file, "%d", &face_id ), "tri id" );
      ref_cell_c2n(ref_cell,3,tri) = face_id;
    }

  ref_cell = ref_grid_qua(ref_grid);
  for( qua = 0; qua < nqua ; qua++ )
    {
      RES( 1, fscanf( file, "%d", &face_id ), "qua id" );
      ref_cell_c2n(ref_cell,4,qua) = face_id;
    }

  ref_cell = ref_grid_tet(ref_grid);
  for( cell = 0; cell < ntet ; cell++ )
    {
      for ( node = 0 ; node < 4 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "tet" );
      nodes[0]--; nodes[1]--; nodes[2]--; nodes[3]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "new tet");
      RES( cell, new_cell, "tet index");
    }

  ref_cell = ref_grid_pyr(ref_grid);
  for( cell = 0; cell < npyr ; cell++ )
    {
      for ( node = 0 ; node < 5 ; node++ )  
	RES( 1, fscanf( file, "%d", &(nodes[node]) ), "pyr" );
      nodes[0]--; nodes[1]--; nodes[2]--; nodes[3]--; nodes[4]--;
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "new pyr");
      RES( cell, new_cell, "pyr index");
    }

  ref_cell = ref_grid_pri(ref_grid);
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

  ref_cell = ref_grid_hex(ref_grid);
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

REF_STATUS ref_import_b8_ugrid( REF_GRID *ref_grid_ptr, char *filename )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_INT node, new_node;
  REF_DBL swapped_dbl;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tri, qua, new_tri, new_qua;
  REF_INT face_id;
  REF_INT node_per, cell, new_cell;

  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  RES( 1, fread( &nnode, sizeof(REF_INT), 1, file ), "nnode" );
  RES( 1, fread( &ntri, sizeof(REF_INT), 1, file ), "ntri" );
  RES( 1, fread( &nqua, sizeof(REF_INT), 1, file ), "nqua" );
  RES( 1, fread( &ntet, sizeof(REF_INT), 1, file ), "ntet" );
  RES( 1, fread( &npyr, sizeof(REF_INT), 1, file ), "npyr" );
  RES( 1, fread( &npri, sizeof(REF_INT), 1, file ), "npri" );
  RES( 1, fread( &nhex, sizeof(REF_INT), 1, file ), "nhex" );

  SWAP_INT(nnode);
  SWAP_INT(ntri);
  SWAP_INT(nqua);
  SWAP_INT(ntet);
  SWAP_INT(npyr);
  SWAP_INT(npri);
  SWAP_INT(nhex);

  for( node=0; node<nnode ; node++ ) 
    {
      RSS( ref_node_add(ref_node, node, &new_node ), "new_node");
      RES( node, new_node, "node index");
      RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "x" );
      SWAP_DBL(swapped_dbl);
      ref_node_xyz( ref_node, 0, new_node ) = swapped_dbl;
      RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "y" );
      SWAP_DBL(swapped_dbl);
      ref_node_xyz( ref_node, 1, new_node ) = swapped_dbl;
      RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "z" );
      SWAP_DBL(swapped_dbl);
      ref_node_xyz( ref_node, 2, new_node ) = swapped_dbl;
    }

  RSS( ref_node_initialize_n_global( ref_node, nnode ), "init glob");

  ref_cell = ref_grid_tri(ref_grid);
  nodes[3] = REF_EMPTY;
  for( tri = 0; tri < ntri ; tri++ )
    {
      node_per = ref_cell_node_per(ref_cell);
      for ( node = 0 ; node < node_per ; node++ )
	{  
	  RES(1, fread( &(nodes[node]), sizeof(REF_INT), 1, file ), "tri" );
	  SWAP_INT(nodes[node]);
	  nodes[node]--;
	}
      RSS( ref_cell_add(ref_cell, nodes, &new_tri ), "new tri");
      RES( tri, new_tri, "tri index");
    }

  ref_cell = ref_grid_qua(ref_grid);
  nodes[4] = REF_EMPTY;
  for( qua = 0; qua < nqua ; qua++ )
    {
      node_per = ref_cell_node_per(ref_cell);
      for ( node = 0 ; node < node_per ; node++ )
	{  
	  RES(1, fread( &(nodes[node]), sizeof(REF_INT), 1, file ), "qua" );
	  SWAP_INT(nodes[node]);
	  nodes[node]--;
	}
      RSS( ref_cell_add(ref_cell, nodes, &new_qua ), "new qua");
      RES( qua, new_qua, "qua index");
    }

  ref_cell = ref_grid_tri(ref_grid);
  for( tri = 0; tri < ntri ; tri++ )
    {
      RES(1, fread( &face_id, sizeof(REF_INT), 1, file ), "tri" );
      SWAP_INT(face_id);
      ref_cell_c2n(ref_cell,3,tri) = face_id;
    }

  ref_cell = ref_grid_qua(ref_grid);
  for( qua = 0; qua < nqua ; qua++ )
    {
      RES(1, fread( &face_id, sizeof(REF_INT), 1, file ), "qua" );
      SWAP_INT(face_id);
      ref_cell_c2n(ref_cell,4,qua) = face_id;
    }

  ref_cell = ref_grid_tet(ref_grid);
  for( cell = 0; cell < ntet ; cell++ )
    {
      node_per = ref_cell_node_per(ref_cell);
      for ( node = 0 ; node < node_per ; node++ )  
	{
	  RES(1, fread( &(nodes[node]), sizeof(REF_INT), 1, file ), "tet" );
	  SWAP_INT(nodes[node]);
	  nodes[node]--;
	}
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "tet cell");
      RES( cell, new_cell, "tet index");
    }

  ref_cell = ref_grid_pyr(ref_grid);
  for( cell = 0; cell < npyr ; cell++ )
    {
      node_per = ref_cell_node_per(ref_cell);
      for ( node = 0 ; node < node_per ; node++ )  
	{
	  RES(1, fread( &(nodes[node]), sizeof(REF_INT), 1, file ), "pyr" );
	  SWAP_INT(nodes[node]);
	  nodes[node]--;
	}
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "pyr cell");
      RES( cell, new_cell, "pyr index");
    }

  ref_cell = ref_grid_pri(ref_grid);
  for( cell = 0; cell < npri ; cell++ )
    {
      node_per = ref_cell_node_per(ref_cell);
      for ( node = 0 ; node < node_per ; node++ )  
	{
	  RES(1, fread( &(nodes[node]), sizeof(REF_INT), 1, file ), "pri" );
	  SWAP_INT(nodes[node]);
	  nodes[node]--;
	}
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "pri cell");
      RES( cell, new_cell, "pri index");
    }

  ref_cell = ref_grid_hex(ref_grid);
  for( cell = 0; cell < nhex ; cell++ )
    {
      node_per = ref_cell_node_per(ref_cell);
      for ( node = 0 ; node < node_per ; node++ )  
	{
	  RES(1, fread( &(nodes[node]), sizeof(REF_INT), 1, file ), "hex" );
	  SWAP_INT(nodes[node]);
	  nodes[node]--;
	}
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "hex cell");
      RES( cell, new_cell, "hex index");
    }

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_import_msh( REF_GRID *ref_grid_ptr, char *filename )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  FILE *file;
  char line[1024];
  REF_INT dummy;
  REF_DBL x,  y;
  REF_INT nnode, node, new_node;
  REF_INT nedge, edge, n0, n1, n2, id;
  REF_INT ntri, tri;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER], new_cell;

  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );
  
  while (!feof(file))
    {
      REIS( 1, fscanf( file, "%s", line), "line read failed");

      if ( 0 == strcmp("Vertices",line))
	{
	  REIS( 1, fscanf(file, "%d", &nnode), "read nnode" );
	  for (node=0;node<nnode;node++)
	    {
	      REIS( 3, fscanf(file, "%lf %lf %d",&x, &y, &dummy ), "read xy" );
	      RSS( ref_node_add( ref_node, node, &new_node ), "add node");
	      ref_node_xyz(ref_node,0,new_node) = x;
	      ref_node_xyz(ref_node,1,new_node) = 0.0;
	      ref_node_xyz(ref_node,2,new_node) = y;
	    }
	  for (node=0;node<nnode;node++)
	    {
	      RSS( ref_node_add( ref_node, nnode+node, &new_node ), "add node");
	      ref_node_xyz(ref_node,0,new_node) = 
		ref_node_xyz(ref_node,0,node);
	      ref_node_xyz(ref_node,1,new_node) = 1.0;
	      ref_node_xyz(ref_node,2,new_node) = 
		ref_node_xyz(ref_node,2,node);
	    }
	  RSS(ref_node_initialize_n_global( ref_node, 2 * nnode ), "init glob");
	}

      if ( 0 == strcmp("Edges",line))
	{
	  REIS( 1, fscanf(file, "%d", &nedge), "read nedge" );
	  for (edge=0;edge<nedge;edge++)
	    {
	      REIS( 3, fscanf(file, "%d %d %d",&n0, &n1, &id ), "read edge" );
	      n0--; n1--;
	      nodes[0]=n0;
	      nodes[1]=n1;
	      nodes[2]=n1+nnode;
	      nodes[3]=n0+nnode;
	      nodes[4]=id;
	      RSS( ref_cell_add( ref_grid_qua(ref_grid), nodes, &new_cell ), 
		   "quad face for an edge");
	    }
	}

      if ( 0 == strcmp("Triangles",line))
	{
	  REIS( 1, fscanf(file, "%d", &ntri), "read ntri" );
	  for (tri=0;tri<ntri;tri++)
	    {
	      REIS( 4, fscanf(file, "%d %d %d %d",&n0, &n1, &n2, &dummy ), 
		    "read tri" );
	      n0--; n1--; n2--;
	      nodes[0]=n0+nnode;
	      nodes[1]=n1+nnode;
	      nodes[2]=n2+nnode;
	      nodes[3]=100;
	      RSS( ref_cell_add( ref_grid_tri(ref_grid), nodes, &new_cell ), 
		   "tri face for tri");
	      nodes[0]=n0;
	      nodes[1]=n2;
	      nodes[2]=n1;
	      nodes[3]=101;
	      RSS( ref_cell_add( ref_grid_tri(ref_grid), nodes, &new_cell ), 
		   "tri face for tri");
	      nodes[0]=n0+nnode;
	      nodes[1]=n1+nnode;
	      nodes[2]=n2+nnode;
	      nodes[3]=n0;
	      nodes[4]=n1;
	      nodes[5]=n2;
	      RSS( ref_cell_add( ref_grid_pri(ref_grid), nodes, &new_cell ), 
		   "prism for tri");
	    }
	  return REF_SUCCESS;
	}

    } 

  return REF_IMPLEMENT;
}

REF_STATUS ref_import_mapbc( REF_DICT *ref_dict_ptr, char *filename )
{
  FILE *file;
  REF_INT n, i;
  REF_INT key, value;
  REF_DICT ref_dict;

  RSS( ref_dict_create( ref_dict_ptr ), "create grid");
  ref_dict = *ref_dict_ptr;

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  RES( 1, fscanf( file, "%d", &n ), "n" );
  for (i=0;i<n;i++)
    {
      RES( 1, fscanf( file, "%d", &key ),   "read mapbc line faceid" );
      RES( 1, fscanf( file, "%d%*[^\n]", &value ), "read mapbc line bc type" );
      RSS( ref_dict_store( ref_dict, key, value ), "store" );
    }

  return REF_SUCCESS;
}
