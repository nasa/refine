
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

  RSS( ref_node_initialize_n_global( ref_node, nnode ), "init nnodesg");

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

  /* ghost xyz */

  RSS( ref_node_ghost_real( ref_node ), "ghost real");

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
	  
	  RSS( ref_sort_heap_int( section_size, dest, order ), "heap" );

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
					     sent_c2n, sent_part,
					     ref_mpi_id ),"many glob");

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
					   sent_c2n, sent_part,
					   ref_mpi_id ), "many glob");

	    ref_free( sent_part );
	  }
      } while ( elements_to_receive != end_of_message );
    }

  free(sent_c2n);

  RSS( ref_migrate_shufflin_cell( ref_node, ref_cell ), "fill ghosts");

  return REF_SUCCESS;
}

REF_STATUS ref_part_metric( REF_NODE ref_node, char *filename )
{
  FILE *file;
  REF_INT chunk;
  REF_DBL *metric;
  REF_INT nnode_read, section_size;
  REF_INT node, local, global, im;

  file = NULL;
  if ( ref_mpi_master )
    {
      file = fopen(filename,"r");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );
    }

  chunk = MAX(100000, ref_node_n_global(ref_node)/ref_mpi_n);
  chunk = MIN( chunk, ref_node_n_global(ref_node) );

  ref_malloc_init( metric, 6*chunk, REF_DBL, -1.0 );
  
  nnode_read = 0;
  while ( nnode_read < ref_node_n_global(ref_node) )
    {
      section_size = MIN(chunk,ref_node_n_global(ref_node)-nnode_read);
      if ( ref_mpi_master )
	{
	  for (node=0;node<section_size;node++)
	    REIS( 6, fscanf( file, "%lf %lf %lf %lf %lf %lf", 
			     &(metric[0+6*node]),
			     &(metric[1+6*node]),
			     &(metric[2+6*node]),
			     &(metric[3+6*node]),
			     &(metric[4+6*node]),
			     &(metric[5+6*node]) ), "metric read error" );
	  RSS( ref_mpi_bcast( metric, 6*chunk, REF_DBL_TYPE ), "bcast" );
	}
      else
	{
	  RSS( ref_mpi_bcast( metric, 6*chunk, REF_DBL_TYPE ), "bcast" );
	}
      for (node=0;node<section_size;node++)
	{
	  global = node + nnode_read;
	  RXS( ref_node_local( ref_node, global, &local ), 
	       REF_NOT_FOUND, "local" );
	  if ( REF_EMPTY != local )
	    for(im=0;im<6;im++)
	      ref_node_metric(ref_node,im,local) = metric[im+6*node];
	}
      nnode_read += section_size;
    }

  ref_free( metric );

  return REF_SUCCESS;
}

REF_STATUS ref_part_ratio( REF_NODE ref_node, REF_DBL *ratio, char *filename )
{
  FILE *file;
  REF_INT chunk;
  REF_DBL *data;
  REF_INT nnode_read, section_size;
  REF_INT node, local, global;

  file = NULL;
  if ( ref_mpi_master )
    {
      file = fopen(filename,"r");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );
    }

  chunk = MAX(100000, ref_node_n_global(ref_node)/ref_mpi_n);
  chunk = MIN( chunk, ref_node_n_global(ref_node) );

  ref_malloc_init( data, chunk, REF_DBL, -1.0 );
  
  nnode_read = 0;
  while ( nnode_read < ref_node_n_global(ref_node) )
    {
      section_size = MIN(chunk,ref_node_n_global(ref_node)-nnode_read);
      if ( ref_mpi_master )
	{
	  for (node=0;node<section_size;node++)
	    REIS( 1, fscanf( file, "%lf", 
			     &(data[node]) ), "data ratio read error" );
	  RSS( ref_mpi_bcast( data, chunk, REF_DBL_TYPE ), "bcast" );
	}
      else
	{
	  RSS( ref_mpi_bcast( data, chunk, REF_DBL_TYPE ), "bcast" );
	}
      for (node=0;node<section_size;node++)
	{
	  global = node + nnode_read;
	  RXS( ref_node_local( ref_node, global, &local ), 
	       REF_NOT_FOUND, "local" );
	  if ( REF_EMPTY != local )
	    ratio[local] = data[node];
	}
      nnode_read += section_size;
    }

  ref_free( data );

  return REF_SUCCESS;
}

