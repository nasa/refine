
#include <stdlib.h>
#include <stdio.h>

#include "ref_part.h"
#include "ref_mpi.h"
#include "ref_endian.h"
#include "ref_sort.h"

#include "ref_malloc.h"

#include "ref_migrate.h"

#include "ref_export.h"

REF_STATUS ref_part_b8_ugrid( REF_GRID *ref_grid_ptr, char *filename )
{
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_INT node, new_node;
  REF_DBL swapped_dbl;
  REF_INT part;
  REF_INT n;
  REF_DBL *xyz;

  long conn_offset, faceid_offset;

  REF_GRID ref_grid;
  REF_NODE ref_node;

  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  /* header */

  file = NULL;
  if ( ref_mpi_master )
    {
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
    }

  RSS( ref_mpi_bcast( &nnode, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &ntri, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &nqua, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &ntet, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &npyr, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &npri, 1, REF_INT_TYPE ), "bcast" ); 
  RSS( ref_mpi_bcast( &nhex, 1, REF_INT_TYPE ), "bcast" ); 

  /* nodes */

  if ( ref_mpi_master )
    {
      part = 0;
      for (node=0;node<ref_part_first( nnode, ref_mpi_n, 1 ); node++)
	{
	  RSS( ref_node_add(ref_node, node, &new_node ), "new_node");
	  ref_node_part(ref_node,new_node) = ref_mpi_id;
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
      for ( part = 1; part<ref_mpi_n ; part++ )
	{
	  n = ref_part_first( nnode, ref_mpi_n, part+1 )
            - ref_part_first( nnode, ref_mpi_n, part );
	  RSS( ref_mpi_send( &n, 1, REF_INT_TYPE, part ), "send" );
	  if ( n > 0 )
	    {
	      ref_malloc( xyz, 3*n, REF_DBL);
	      for (node=0;node<n; node++)
		{
		  RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "x" );
		  SWAP_DBL(swapped_dbl);
		  xyz[0+3*node] = swapped_dbl;
		  RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "y" );
		  SWAP_DBL(swapped_dbl);
		  xyz[1+3*node] = swapped_dbl;
		  RES(1, fread( &swapped_dbl, sizeof(REF_DBL), 1, file ), "z" );
		  SWAP_DBL(swapped_dbl);
		  xyz[2+3*node] = swapped_dbl;
		}
	      RSS( ref_mpi_send( xyz, 3*n, REF_DBL_TYPE, part ), "send" );
	      free(xyz);
	    }
	}
    }
  else
    {
      RSS( ref_mpi_recv( &n, 1, REF_INT_TYPE, 0 ), "recv" );
      if ( n > 0 )
	{
	  ref_malloc( xyz, 3*n, REF_DBL);
	  RSS( ref_mpi_recv( xyz, 3*n, REF_DBL_TYPE, 0 ), "recv" );
	  for (node=0;node<n; node++)
	    {
	      RSS( ref_node_add(ref_node, 
				node+ref_part_first( nnode, 
						     ref_mpi_n, 
						     ref_mpi_id ),
				&new_node ), "new_node");
	      ref_node_part(ref_node,new_node) = ref_mpi_id;
	      ref_node_xyz( ref_node, 0, new_node ) = xyz[0+3*node];
	      ref_node_xyz( ref_node, 1, new_node ) = xyz[1+3*node];
	      ref_node_xyz( ref_node, 2, new_node ) = xyz[2+3*node];
	    }
	  free(xyz);
	}

    }

  if ( 0 < ntri )
    {
      conn_offset = 4*7
	+ 8*3*nnode;
      faceid_offset = 4*7
	+ 8*3*nnode
	+ 4*3*ntri
	+ 4*4*nqua;
      RSS( ref_part_b8_ugrid_cell( ref_grid_tri(ref_grid), ntri, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "tri" );
    }

  if ( 0 < nqua )
    {
      conn_offset = 4*7
	+ 8*3*nnode
	+ 4*3*ntri;
      faceid_offset = 4*7
	+ 8*3*nnode
	+ 4*4*ntri
	+ 4*4*nqua;
      RSS( ref_part_b8_ugrid_cell( ref_grid_qua(ref_grid), nqua, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "qua" );
    }

  if ( 0 < ntet )
    {
      conn_offset = 4*7
	+ 8*3*nnode
	+ 4*4*ntri
	+ 4*5*nqua;
      faceid_offset = REF_EMPTY;
      RSS( ref_part_b8_ugrid_cell( ref_grid_tet(ref_grid), ntet, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "tet" );
    }

  if ( 0 < npyr )
    {
      conn_offset = 4*7
	+ 8*3*nnode
	+ 4*4*ntri
	+ 4*5*nqua
	+ 4*4*ntet;
      faceid_offset = REF_EMPTY;
      RSS( ref_part_b8_ugrid_cell( ref_grid_pyr(ref_grid), npyr, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "pyr" );
    }

  if ( 0 < npri )
    {
      conn_offset = 4*7
	+ 8*3*nnode
	+ 4*4*ntri
	+ 4*5*nqua
	+ 4*4*ntet
	+ 4*5*npyr;
      faceid_offset = REF_EMPTY;
      RSS( ref_part_b8_ugrid_cell( ref_grid_pri(ref_grid), npri, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "pri" );
    }

  if ( 0 < nhex )
    {
      conn_offset = 4*7
	+ 8*3*nnode
	+ 4*4*ntri
	+ 4*5*nqua
	+ 4*4*ntet
	+ 4*5*npyr
	+ 4*6*npri;
      faceid_offset = REF_EMPTY;
      RSS( ref_part_b8_ugrid_cell( ref_grid_hex(ref_grid), nhex, 
				   ref_node, nnode, 
				   file, conn_offset, faceid_offset ), "hex" );
    }

  if ( ref_mpi_master ) fclose(file);

  ref_node_n_global(ref_node) = nnode;

  /* ghost xyz */

  RSS( ref_part_ghost_xyz( ref_grid ), "ghost xyz");

  return REF_SUCCESS;
}

REF_STATUS ref_part_b8_ugrid_cell( REF_CELL ref_cell, REF_INT ncell,
				   REF_NODE ref_node, REF_INT nnode,
				   FILE *file, 
				   long conn_offset, 
				   long faceid_offset )
{
  REF_INT ncell_read;
  REF_INT chunk;
  REF_INT end_of_message = REF_EMPTY;
  REF_INT elements_to_receive;
  REF_INT *c2n;
  REF_INT *sent_c2n;
  REF_INT *dest;
  REF_INT *order;
  REF_INT *sent_part;
  REF_INT *elements_to_send;
  REF_INT *start_to_send;
  REF_INT node_per, size_per;
  REF_INT section_size;
  REF_INT cell;
  REF_INT part, node;
  REF_INT ncell_keep;

  chunk = MAX(1000000, ncell/ref_mpi_n);

  size_per = ref_cell_size_per(ref_cell);
  node_per = ref_cell_node_per(ref_cell);

  ref_malloc( sent_c2n, size_per*chunk, REF_INT );

  if ( ref_mpi_master )
    {

      ref_malloc( elements_to_send, ref_mpi_n, REF_INT );
      ref_malloc( start_to_send, ref_mpi_n, REF_INT );
      ref_malloc( c2n, size_per*chunk, REF_INT );
      ref_malloc( dest, chunk, REF_INT );
      ref_malloc( order, chunk, REF_INT );

      ncell_read = 0;
      while ( ncell_read < ncell )
	{
	  section_size = MIN(chunk,ncell-ncell_read);

	  fseek(file,conn_offset+(long)(4*node_per*ncell_read),SEEK_SET);
	  for (cell=0;cell<section_size;cell++)
	    {
	      for (node=0;node<node_per;node++)
		{
		  RES(1, fread( &(c2n[node+size_per*cell]), 
				sizeof(REF_INT), 1, file ), "cn" );
		  SWAP_INT(c2n[node+size_per*cell]);
		  c2n[node+size_per*cell]--;
		}
	    }
	  if ( ref_cell_last_node_is_an_id(ref_cell) )
	    {
	      fseek(file,faceid_offset+(long)(4*ncell_read),SEEK_SET);
	      for (cell=0;cell<section_size;cell++)
		{
		  node = node_per;
		  RES(1, fread( &(c2n[node+size_per*cell]), 
				sizeof(REF_INT), 1, file ), "cn" );
		  SWAP_INT(c2n[node+size_per*cell]);
		}
	    }

	  ncell_read += section_size;

	  for (cell=0;cell<section_size;cell++)
	    dest[cell] = ref_part_implicit( nnode, ref_mpi_n, 
					    c2n[size_per*cell] );
	  
	  RSS( ref_sort_heap( section_size, dest, order ), "heap" );

	  for (cell=0;cell<section_size;cell++)
	    {
	      for (node=0;node<size_per;node++)
		sent_c2n[node+size_per*cell] = c2n[node+size_per*order[cell]];
	      dest[cell] = ref_part_implicit( nnode, ref_mpi_n, 
					      sent_c2n[size_per*cell] );
	    }

	  /* master keepers */
	  
	  for ( part = 0; part< ref_mpi_n;part++ )
	    elements_to_send[part] = 0;
	  for ( part = 0; part< ref_mpi_n;part++ )
	    start_to_send[part] = REF_EMPTY;

	  for (cell=section_size-1; cell >= 0; cell--)
	    {
	      elements_to_send[dest[cell]]++;
	      start_to_send[dest[cell]] = cell;
	    }

	  ncell_keep = elements_to_send[0];
	  if ( 0 < ncell_keep )
	    { 
	      ref_malloc_init( sent_part, size_per*ncell_keep, 
			       REF_INT, REF_EMPTY );

	      for (cell=0;cell<ncell_keep;cell++)
		for (node=0;node<node_per;node++)
		  sent_part[node+size_per*cell] =
		    ref_part_implicit( nnode, ref_mpi_n, 
				       sent_c2n[node+size_per*cell] );

	      RSS( ref_cell_add_many_global( ref_cell, ref_node,
					     ncell_keep, 
					     sent_c2n, sent_part ),"many glob");

	      ref_free(sent_part);
	    }

	  for ( part = 1; part< ref_mpi_n;part++ )
	    if ( 0 < elements_to_send[part] ) 
	      {
		RSS( ref_mpi_send( &(elements_to_send[part]), 
				   1, REF_INT_TYPE, part ), "send" );
		RSS( ref_mpi_send( &(sent_c2n[size_per*start_to_send[part]]), 
				   size_per*elements_to_send[part], 
				   REF_INT_TYPE, part ), "send" );
	      }

	}

      ref_free(order);
      ref_free(dest);
      ref_free(c2n);
      ref_free(start_to_send);
      ref_free(elements_to_send);

      /* signal we are done */
      for ( part = 1; part<ref_mpi_n ; part++ )
	RSS( ref_mpi_send( &end_of_message, 1, REF_INT_TYPE, part ), "send" );

    }  
  else
    {
      do {
	RSS( ref_mpi_recv( &elements_to_receive, 1, REF_INT_TYPE, 0 ), "recv" );
	if ( elements_to_receive > 0 )
	  {
	    RSS( ref_mpi_recv( sent_c2n, size_per*elements_to_receive, 
			       REF_INT_TYPE, 0 ), "send" );

	    ref_malloc_init( sent_part, size_per*elements_to_receive, 
			     REF_INT, REF_EMPTY );

	    for (cell=0;cell<elements_to_receive;cell++)
	      for (node=0;node<node_per;node++)
		sent_part[node+size_per*cell] =
		  ref_part_implicit( nnode, ref_mpi_n, 
				     sent_c2n[node+size_per*cell] );

	    RSS( ref_cell_add_many_global( ref_cell, ref_node,
					   elements_to_receive, 
					   sent_c2n, sent_part ), "many glob");

	    ref_free( sent_part );
	  }
      } while ( elements_to_receive != end_of_message );
    }

  free(sent_c2n);

  RSS( ref_migrate_shufflin_cell( ref_node, ref_cell ), "fill ghosts");

  return REF_SUCCESS;
}

REF_STATUS ref_part_ghost_xyz( REF_GRID ref_grid )
{
  REF_NODE ref_node;
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT *a_global, *b_global;
  REF_INT part, node;
  REF_INT *a_next;
  REF_DBL *a_xyz, *b_xyz;
  REF_INT local;

  if ( 1 == ref_mpi_n ) return REF_SUCCESS;

  ref_node = ref_grid_node(ref_grid);

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
	a_next[ref_node_part(ref_node,node)]++;
      }

  RSS( ref_mpi_alltoallv( a_global, a_size, b_global, b_size, 
			  1, REF_INT_TYPE ), 
       "alltoallv global");

  for (node=0;node<b_total;node++)
    {
      RSS( ref_node_local( ref_node, b_global[node], &local ), "g2l");
      b_xyz[0+3*node] = ref_node_xyz(ref_node,0,local);
      b_xyz[1+3*node] = ref_node_xyz(ref_node,1,local);
      b_xyz[2+3*node] = ref_node_xyz(ref_node,2,local);
    }

  RSS( ref_mpi_alltoallv( b_xyz, b_size, a_xyz, a_size, 
			  3, REF_DBL_TYPE ), 
       "alltoallv global");

  for (node=0;node<a_total;node++)
    {
      RSS( ref_node_local( ref_node, a_global[node], &local ), "g2l");
      ref_node_xyz(ref_node,0,local) = a_xyz[0+3*node];
      ref_node_xyz(ref_node,1,local) = a_xyz[1+3*node];
      ref_node_xyz(ref_node,2,local) = a_xyz[2+3*node];
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

REF_STATUS ref_part_ghost_int( REF_GRID ref_grid, REF_INT *scalar )
{
  REF_NODE ref_node;
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT *a_global, *b_global;
  REF_INT part, node;
  REF_INT *a_next;
  REF_INT *a_scalar, *b_scalar;
  REF_INT local;

  if ( 1 == ref_mpi_n ) return REF_SUCCESS;

  ref_node = ref_grid_node(ref_grid);

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
  ref_malloc( a_scalar, a_total, REF_INT );

  b_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    b_total += b_size[part];
  ref_malloc( b_global, b_total, REF_INT );
  ref_malloc( b_scalar, b_total, REF_INT );

  ref_malloc( a_next, ref_mpi_n, REF_INT );
  a_next[0] = 0;
  for ( part = 1; part<ref_mpi_n ; part++ )
    a_next[part] = a_next[part-1]+a_size[part-1];

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      {
	part = ref_node_part(ref_node,node);
	a_global[a_next[part]] = ref_node_global(ref_node,node);
	a_next[ref_node_part(ref_node,node)]++;
      }

  RSS( ref_mpi_alltoallv( a_global, a_size, b_global, b_size, 
			  1, REF_INT_TYPE ), 
       "alltoallv global");

  for (node=0;node<b_total;node++)
    {
      RSS( ref_node_local( ref_node, b_global[node], &local ), "g2l");
      b_scalar[node] = scalar[local];
    }

  RSS( ref_mpi_alltoallv( b_scalar, b_size, a_scalar, a_size, 
			  1, REF_INT_TYPE ), 
       "alltoallv global");

  for (node=0;node<a_total;node++)
    {
      RSS( ref_node_local( ref_node, a_global[node], &local ), "g2l");
      scalar[local] = a_scalar[node];
    }

  free(a_next);
  free(b_scalar);
  free(b_global);
  free(a_scalar);
  free(a_global);
  free(b_size);
  free(a_size);

  return REF_SUCCESS;  
}
