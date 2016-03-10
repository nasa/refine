
#include <stdlib.h>
#include <stdio.h>

#include <string.h>

#include "ref_import.h"

#include "ref_malloc.h"
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

  for(i=0;i<8;i++)
    {
      RES( 1, fread( &i4, sizeof(int), 1, file ), "int" );
      i4_swapped = i4; SWAP_INT(i4_swapped);
      printf(" %d: %d (%d swapped) ints\n",i,i4,i4_swapped);
    }

  printf(" --\n");

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

  if( strcmp(&filename[end_of_string-10],".lb8.ugrid") == 0 ) 
    {
      RSS( ref_import_lb8_ugrid( ref_grid_ptr, filename ), "lb8_ugrid failed");
    } 
  else 
    if( strcmp(&filename[end_of_string-9],".b8.ugrid") == 0 ) 
      {
	RSS( ref_import_b8_ugrid( ref_grid_ptr, filename ), "b8_ugrid failed");
      } 
    else 
      if( strcmp(&filename[end_of_string-9],".r8.ugrid") == 0 ) 
	{
	  RSS( ref_import_r8_ugrid( ref_grid_ptr, filename ), "r8_ugrid failed");
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
	      if( strcmp(&filename[end_of_string-6],".meshb") == 0 ) 
		{
		  RSS( ref_import_meshb( ref_grid_ptr, filename ), 
		       "meshb failed");
		} 
	      else 
		{
		  printf("%s: %d: %s %s\n",__FILE__,__LINE__,
			 "input file name extension unknown", filename);
		  RSS( REF_FAILURE, "unknown file extension");
		}
  
  ref_grid_guess_twod_status( *ref_grid_ptr );

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

static REF_STATUS ref_import_bin_ugrid_bound_c2n( REF_CELL ref_cell,
						  REF_INT ncell,
						  FILE *file,
						  REF_BOOL swap )
{
  REF_INT node_per, max_chunk, nread, chunk, cell, node, new_cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT *c2n;

  if ( 0 < ncell )
    {
      node_per = ref_cell_node_per(ref_cell);
      max_chunk = MIN(1000000, ncell);
      ref_malloc( c2n, node_per*max_chunk, REF_INT);
      nread = 0;
      while ( nread < ncell )
	{
	  chunk = MIN(max_chunk, ncell-nread);
	  REIS((size_t)(node_per*chunk),
	       fread( c2n, sizeof(REF_INT), node_per*chunk, file ), "c2n" );
	  for( cell=0; cell<chunk ; cell++ )
	    {
	      nodes[node_per] = REF_EMPTY;
	      for ( node = 0 ; node < node_per ; node++ )
		{
		  nodes[node] = c2n[node+node_per*cell];
		  if (swap) SWAP_INT(nodes[node]);
		  nodes[node]--;
		}
	      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "new cell");
	      RES( cell+nread, new_cell, "cell index");
	    }
	  nread += ncell;
	}
      ref_free( c2n );
    }

  return REF_SUCCESS;
}

static REF_STATUS ref_import_bin_ugrid_bound_tag( REF_CELL ref_cell,
						  REF_INT ncell,
						  FILE *file,
						  REF_BOOL swap )
{
  REF_INT node_per, max_chunk, nread, chunk, cell;
  REF_INT *tag;

  if ( 0 < ncell )
    {
      node_per = ref_cell_node_per(ref_cell);
      max_chunk = MIN(1000000, ncell);
      ref_malloc( tag, max_chunk, REF_INT);
      nread = 0;
      while ( nread < ncell )
	{
	  chunk = MIN(max_chunk, ncell-nread);
	  REIS((size_t)(chunk),
	       fread( tag, sizeof(REF_INT), chunk, file ), "tags" );
	  for( cell=0; cell<chunk ; cell++ )
	    {
	      if (swap) SWAP_INT(tag[cell]);
	      ref_cell_c2n(ref_cell,node_per,cell) = tag[cell];
	    }
	  nread += ncell;
	}
      ref_free( tag );
    }

  return REF_SUCCESS;
}

static REF_STATUS ref_import_bin_ugrid( REF_GRID *ref_grid_ptr, char *filename,
					REF_BOOL swap )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_INT node, new_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell, new_cell;

  REF_INT max_chunk, nread, chunk, ixyz;
  REF_DBL *xyz;

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

  if (swap) SWAP_INT(nnode);
  if (swap) SWAP_INT(ntri);
  if (swap) SWAP_INT(nqua);
  if (swap) SWAP_INT(ntet);
  if (swap) SWAP_INT(npyr);
  if (swap) SWAP_INT(npri);
  if (swap) SWAP_INT(nhex);

  /* large block reads reccomended for IO performance */
  max_chunk = MIN(1000000, nnode);
  ref_malloc( xyz, 3*max_chunk, REF_DBL);
  nread = 0;
  while ( nread < nnode )
    {
      chunk = MIN(max_chunk, nnode-nread);
      REIS((size_t)(3*chunk),
	   fread( xyz, sizeof(REF_DBL), 3*chunk, file ), "xyz" );
      if (swap)
	for( ixyz=0; ixyz<3*chunk ; ixyz++ )
	  SWAP_DBL(xyz[ixyz]);
      for( node=0; node<chunk ; node++ )
	{
	  RSS( ref_node_add(ref_node, node+nread, &new_node ), "new_node");
	  ref_node_xyz( ref_node, 0, new_node ) = xyz[0+3*node];
	  ref_node_xyz( ref_node, 1, new_node ) = xyz[1+3*node];
	  ref_node_xyz( ref_node, 2, new_node ) = xyz[2+3*node];
	}
      nread += chunk;
    }
  ref_free( xyz );

  RSS( ref_node_initialize_n_global( ref_node, nnode ), "init glob");

  RSS( ref_import_bin_ugrid_bound_c2n( ref_grid_tri(ref_grid), ntri,
				       file, swap ), "tri face nodes");
  RSS( ref_import_bin_ugrid_bound_c2n( ref_grid_qua(ref_grid), nqua,
				       file, swap ), "qua face nodes");

  RSS( ref_import_bin_ugrid_bound_tag( ref_grid_tri(ref_grid), ntri,
				       file, swap ), "tri face tags");
  RSS( ref_import_bin_ugrid_bound_tag( ref_grid_qua(ref_grid), nqua,
				       file, swap ), "tri face tags");

  ref_cell = ref_grid_tet(ref_grid);
  for( cell = 0; cell < ntet ; cell++ )
    {
      node_per = ref_cell_node_per(ref_cell);
      for ( node = 0 ; node < node_per ; node++ )  
	{
	  RES(1, fread( &(nodes[node]), sizeof(REF_INT), 1, file ), "tet" );
	  if (swap) SWAP_INT(nodes[node]);
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
	  if (swap) SWAP_INT(nodes[node]);
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
	  if (swap) SWAP_INT(nodes[node]);
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
	  if (swap) SWAP_INT(nodes[node]);
	  nodes[node]--;
	}
      RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "hex cell");
      RES( cell, new_cell, "hex index");
    }

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_import_lb8_ugrid( REF_GRID *ref_grid_ptr, char *filename )
{
  RSS( ref_import_bin_ugrid( ref_grid_ptr, filename, REF_FALSE ),
       "import bin ugrid unswapped");
  return REF_SUCCESS;
}

REF_STATUS ref_import_b8_ugrid( REF_GRID *ref_grid_ptr, char *filename )
{
  RSS( ref_import_bin_ugrid( ref_grid_ptr, filename, REF_TRUE ),
       "import bin ugrid swapped");
  return REF_SUCCESS;
}

REF_STATUS ref_import_r8_ugrid( REF_GRID *ref_grid_ptr, char *filename )
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
  REF_INT fortran_record_size;

  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "nnode" );
  SWAP_INT(fortran_record_size);
  REIS( 7*4, fortran_record_size, "header start record size" );

  RES( 1, fread( &nnode, sizeof(REF_INT), 1, file ), "nnode" );
  RES( 1, fread( &ntri, sizeof(REF_INT), 1, file ), "ntri" );
  RES( 1, fread( &nqua, sizeof(REF_INT), 1, file ), "nqua" );
  RES( 1, fread( &ntet, sizeof(REF_INT), 1, file ), "ntet" );
  RES( 1, fread( &npyr, sizeof(REF_INT), 1, file ), "npyr" );
  RES( 1, fread( &npri, sizeof(REF_INT), 1, file ), "npri" );
  RES( 1, fread( &nhex, sizeof(REF_INT), 1, file ), "nhex" );

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "nnode" );
  SWAP_INT(fortran_record_size);
  REIS( 7*4, fortran_record_size, "header end record size" );

  SWAP_INT(nnode);
  SWAP_INT(ntri);
  SWAP_INT(nqua);
  SWAP_INT(ntet);
  SWAP_INT(npyr);
  SWAP_INT(npri);
  SWAP_INT(nhex);

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "nnode" );
  SWAP_INT(fortran_record_size);
  REIS( nnode*3*8 +
	ntri*4*4 +
	nqua*5*4 +
	ntet*4*4 +
	npyr*5*4 +
	npri*6*4 +
	nhex*8*4 , 
	fortran_record_size, "block start record size" );

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

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "nnode" );
  SWAP_INT(fortran_record_size);
  REIS( nnode*3*8 +
	ntri*4*4 +
	nqua*5*4 +
	ntet*4*4 +
	npyr*5*4 +
	npri*6*4 +
	nhex*8*4 , 
	fortran_record_size, "block end record size" );

  fclose(file);

  return REF_SUCCESS;
}

/*
                           2 
                        ./   \
                     .  / nt2 \
                  .    /       \
               .      /nt0  nt1 \
            .        0-ne0---ne1-1  y=1 second plane, face 1
         5       .           .
       /   \  .           .
      /  t1 \          .
     /   .   \      .
    / t0   t2 \  .
   3--e0----e1-4   y=0 first plane, face 2

 */

REF_STATUS ref_import_msh( REF_GRID *ref_grid_ptr, char *filename )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  FILE *file;
  char line[1024];
  REF_INT dummy, row;
  REF_DBL x, y, z;
  REF_INT nnode, node, new_node;
  REF_INT nedge, edge, n0, n1, n2, id;
  REF_INT ntri, tri;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER], new_cell;
  REF_INT status;
  REF_INT elem, nelem, type, flag, three, zero;

  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );
  
  while (!feof(file))
    {
      status = fscanf( file, "%s", line);
      if ( EOF == status ) return REF_SUCCESS;
      REIS( 1, status, "line read failed");

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
	      nodes[3]=1;
	      RSS( ref_cell_add( ref_grid_tri(ref_grid), nodes, &new_cell ), 
		   "tri face for tri");
	      nodes[0]=n0;
	      nodes[1]=n2;
	      nodes[2]=n1;
	      nodes[3]=2;
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
	}

      if ( 0 == strcmp("$Nodes",line))
	{
	  REIS( 1, fscanf(file, "%d", &nnode), "read nnode" );
	  printf("$Nodes\n%d\n",nnode);
	  for (node=0;node<nnode;node++)
	    {
	      REIS( 4, fscanf(file, "%d %lf %lf %lf",&row, &x, &y, &z ), 
		    "read $Nodes xyz" );
	      REIS(node+1,row,"row index missmatch in $Nodes");
	      RSS( ref_node_add( ref_node, node, &new_node ), "add node");
	      ref_node_xyz(ref_node,0,new_node) = x;
	      ref_node_xyz(ref_node,1,new_node) = y;
	      ref_node_xyz(ref_node,2,new_node) = z;
	    }
	  RSS(ref_node_initialize_n_global( ref_node, nnode ), "init glob");
	}

      if ( 0 == strcmp("$Elements",line))
	{
	  REIS( 1, fscanf(file, "%d", &nelem), "read nelements" );
	  printf("$Elements\n%d\n",nelem);
	  for (elem=0;elem<nelem;elem++)
	    {
	      REIS( 6, fscanf(file, "%d %d %d %d %d %d",
			      &row, &type, &three, &id, &flag, &zero ), 
		    "$Elements description" );
	      REIS(elem+1,row,"row index missmatch in $Elements");
	      switch( type )
		{
		case 5:
		  ref_cell = ref_grid_hex(ref_grid);
		  break;
		case 3:
		  ref_cell = ref_grid_qua(ref_grid);
		  break;
		default:
		   printf("type = %d\n",type);
		   THROW("unknown $Elements type");
		}
	      for (node=0;node<ref_cell_node_per(ref_cell);node++)
		{
		  REIS( 1, fscanf(file, "%d", &(nodes[node]) ),
			"$Elements node" );
		  (nodes[node])--;
		}
	      if ( 3 == type ) nodes[4] = id;
	      RSS( ref_cell_add( ref_cell, nodes, &new_cell ), "add $Element");
	    }
	  ref_grid_inspect(ref_grid);
	}

    } 

  fclose(file);

  return REF_IMPLEMENT;
}

static REF_STATUS meshb_real( FILE *file, REF_INT version, REF_DBL *real )
{
  float temp_float;
  double temp_double;

  if ( 1 == version )
    {
      REIS( 1, fread(&temp_float,sizeof(temp_float), 1, file ), "read float" );
      *real = temp_float;
    }
  else
    {
      REIS( 1, fread(&temp_double,sizeof(temp_double), 1, file),"read double");
      *real = temp_double;
    }

  return REF_SUCCESS;
} 

static REF_STATUS meshb_pos( FILE *file, REF_INT version, REF_INT *pos )
{
  int temp_int;
  long temp_long;

  if ( 3 == version )
    {
      REIS( 1, fread(&temp_long,sizeof(temp_long), 1, file ), "read long" );
      *pos = temp_long;
    }
  else
    {
      REIS( 1, fread(&temp_int,sizeof(temp_int), 1, file),"read double");
      *pos = temp_int;
    }

  return REF_SUCCESS;
} 

REF_STATUS ref_import_meshb( REF_GRID *ref_grid_ptr, char *filename )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  FILE *file;
  REF_INT code, version, dim;
  REF_INT keyword_code, position, next_position, end_position;
  REF_DICT ref_dict;
  REF_INT vertex_keyword, triangle_keyword, edge_keyword, tet_keyword;
  REF_INT nnode, node, new_node;
  REF_INT ntri, tri, nedge, edge, ntet, tet;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER], new_cell;
  REF_INT n0, n1, n2, n3, id;
  REF_BOOL verbose = REF_FALSE;
  
  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  RSS( ref_dict_create( &ref_dict ), "create dict" );

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );
  
  REIS(1, fread((unsigned char *)&code, 4, 1, file), "code");
  REIS(1, code, "code");
  REIS(1, fread((unsigned char *)&version, 4, 1, file), "version");
  if ( version < 1 || 3 < version )
    {
      printf("version %d not supported\n",version);
      THROW("version");
    }
  if (verbose) printf("meshb version %d\n",version);

  position = ftell(file);
  REIS(1, fread((unsigned char *)&keyword_code, 4, 1, file), "keyword code");
  REIS(3, keyword_code, "keyword code");
  RSS( ref_dict_store( ref_dict, keyword_code, position ), "store pos");
  RSS( meshb_pos( file, version, &next_position), "pos");

  REIS(1, fread((unsigned char *)&dim, 4, 1, file), "dim");
  if ( dim < 2 || 3 < dim )
    {
      printf("dim %d not supported\n",dim);
      THROW("dim");
    }
  if (verbose) printf("meshb dim %d\n",dim);

  fseek(file, 0, SEEK_END);
  end_position = ftell(file);

  while ( next_position <= end_position && 0 != next_position )
    {
      position = next_position;
      fseek(file, position, SEEK_SET);
      REIS(1, fread((unsigned char *)&keyword_code, 4, 1, file), 
	   "keyword code");
      RSS( ref_dict_store( ref_dict, keyword_code, position ), "store pos");
      RSS( meshb_pos( file, version, &next_position), "pos");
    }  

  if (verbose) ref_dict_inspect(ref_dict);

  vertex_keyword = 4;
  RSS( ref_dict_value( ref_dict, vertex_keyword, &position), "kw pos");
  fseek(file, (long)position, SEEK_SET);
  REIS(1, fread((unsigned char *)&keyword_code, 4, 1, file), "keyword code");
  REIS(vertex_keyword, keyword_code, "keyword code");
  RSS( meshb_pos( file, version, &next_position), "pos");
  REIS(1, fread((unsigned char *)&nnode, 4, 1, file), "keyword code");
  if (verbose) printf("nnode %d\n",nnode);

  for (node=0;node<nnode;node++)
    {
      RSS( ref_node_add( ref_node, node, &new_node ), "add node");
      if ( 2 == dim )
	{
	  RSS(meshb_real(file, version, 
			 &(ref_node_xyz(ref_node,0,new_node)) ),"x");
	  ref_node_xyz(ref_node,1,new_node) = 0.0;
	  RSS(meshb_real(file, version, 
			 &(ref_node_xyz(ref_node,2,new_node)) ),"y");
	}
      else
	{
	  RSS(meshb_real(file, version, 
			 &(ref_node_xyz(ref_node,0,new_node)) ),"x");
	  RSS(meshb_real(file, version, 
			 &(ref_node_xyz(ref_node,1,new_node)) ),"y");
	  RSS(meshb_real(file, version, 
			 &(ref_node_xyz(ref_node,2,new_node)) ),"z");
	}
      REIS( 1, fread(&(id),sizeof(id), 1, file ), "id" );
      /* ref_node_location(ref_node, node ); */
    }
  REIS( next_position, ftell(file), "end location" );
  if ( 2 == dim )
    for (node=0;node<nnode;node++)
      {
	RSS( ref_node_add( ref_node, nnode+node, &new_node ), "add node");
	ref_node_xyz(ref_node,0,new_node) = 
	  ref_node_xyz(ref_node,0,node);
	ref_node_xyz(ref_node,1,new_node) = 1.0;
	ref_node_xyz(ref_node,2,new_node) = 
	  ref_node_xyz(ref_node,2,node);
      }

  RSS( ref_node_initialize_n_global( ref_node, nnode ), "init glob");

  if ( 2 == dim )
    {
      edge_keyword = 5;
      RSS( ref_dict_value( ref_dict, edge_keyword, &position), "kw pos");
      fseek(file, (long)position, SEEK_SET);
      REIS(1, fread((unsigned char *)&keyword_code, 4, 1, file), 
	   "keyword code");
      REIS(edge_keyword, keyword_code, "keyword code");
      RSS( meshb_pos( file, version, &next_position), "pos");
      REIS(1, fread((unsigned char *)&nedge, 4, 1, file), "keyword code");
      if (verbose) printf("nedge %d\n",nedge);

      for (edge=0;edge<nedge;edge++)
	{
	  REIS( 1, fread(&(n0),sizeof(n0), 1, file ), "n0" );
	  REIS( 1, fread(&(n1),sizeof(n1), 1, file ), "n1" );
	  REIS( 1, fread(&(id),sizeof(id), 1, file ), "id" );
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

  triangle_keyword = 6;
  RSS( ref_dict_value( ref_dict, triangle_keyword, &position), "kw pos");
  fseek(file, (long)position, SEEK_SET);
  REIS(1, fread((unsigned char *)&keyword_code, 4, 1, file), "keyword code");
  REIS(triangle_keyword, keyword_code, "keyword code");
  RSS( meshb_pos( file, version, &next_position), "pos");
  REIS(1, fread((unsigned char *)&ntri, 4, 1, file), "keyword code");
  if (verbose) printf("ntri %d\n",ntri);

  for (tri=0;tri<ntri;tri++)
    {
      REIS( 1, fread(&(n0),sizeof(n0), 1, file ), "n0" );
      REIS( 1, fread(&(n1),sizeof(n1), 1, file ), "n1" );
      REIS( 1, fread(&(n2),sizeof(n2), 1, file ), "n2" );
      REIS( 1, fread(&(id),sizeof(id), 1, file ), "id" );
      n0--; n1--; n2--;
      if ( 2 == dim )
	{
	  nodes[0]=n0+nnode;
	  nodes[1]=n1+nnode;
	  nodes[2]=n2+nnode;
	  nodes[3]=1;
	  RSS( ref_cell_add( ref_grid_tri(ref_grid), nodes, &new_cell ), 
	       "tri face for tri");
	  nodes[0]=n0;
	  nodes[1]=n2;
	  nodes[2]=n1;
	  nodes[3]=2;
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
      else
	{
	  nodes[0]=n0;
	  nodes[1]=n1;
	  nodes[2]=n2;
	  nodes[3]=id;
	  RSS( ref_cell_add( ref_grid_tri(ref_grid), nodes, &new_cell ), 
	       "tri face for tri");
	}
    }
  REIS( next_position, ftell(file), "end location" );

  tet_keyword = 8;
  code = ref_dict_value( ref_dict, tet_keyword, &position);
  RXS( code, REF_NOT_FOUND, "kw pos");
  if ( 3==dim && code != REF_NOT_FOUND )
    {
      fseek(file, (long)position, SEEK_SET);
      REIS(1, fread((unsigned char *)&keyword_code, 4, 1, file), 
	   "keyword code");
      REIS(tet_keyword, keyword_code, "keyword code");
      RSS( meshb_pos( file, version, &next_position), "pos");
      REIS(1, fread((unsigned char *)&ntet, 4, 1, file), "keyword code");
      if (verbose) printf("ntet %d\n",ntet);

      for (tet=0;tet<ntet;tet++)
	{
	  REIS( 1, fread(&(n0),sizeof(n0), 1, file ), "n0" );
	  REIS( 1, fread(&(n1),sizeof(n1), 1, file ), "n1" );
	  REIS( 1, fread(&(n2),sizeof(n2), 1, file ), "n2" );
	  REIS( 1, fread(&(n3),sizeof(n3), 1, file ), "n3" );
	  REIS( 1, fread(&(id),sizeof(id), 1, file ), "id" );
	  n0--; n1--; n2--; n3--;
	  nodes[0]=n0;
	  nodes[1]=n1;
	  nodes[2]=n2;
	  nodes[3]=n3;
	  RSS( ref_cell_add( ref_grid_tet(ref_grid), nodes, &new_cell ), 
	       "tet");
	}
    }
  REIS( next_position, ftell(file), "end location" );

  ref_dict_free( ref_dict );

  fclose(file);

  return REF_SUCCESS;
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

  fclose(file);

  return REF_SUCCESS;
}
