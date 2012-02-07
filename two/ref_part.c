
#include <stdlib.h>
#include <stdio.h>

#include "ref_part.h"
#include "ref_mpi.h"
#include "ref_endian.h"
#include "ref_sort.h"

#include "ref_export.h"

REF_STATUS ref_part_b8_ugrid( REF_GRID *ref_grid_ptr, char *filename )
{
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_INT cell;
  REF_INT node, new_node;
  REF_DBL swapped_dbl;
  REF_INT part;
  REF_INT n;
  REF_DBL *xyz;
  REF_INT end_of_message = REF_EMPTY;
  REF_INT elements_to_receive;
  REF_INT *c2n;
  REF_INT *sent_c2n;
  REF_INT *elements_to_send;
  REF_INT node_per;
  REF_INT ncell, section_size;
  REF_INT all_procs[REF_CELL_MAX_NODE_PER];
  REF_INT unique_procs[REF_CELL_MAX_NODE_PER];
  REF_INT cell_procs;

  REF_INT chunk;

  REF_BOOL needcell;

  REF_INT send_size, new_cell;

  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;

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
	      xyz=(REF_DBL*)malloc(3*n*sizeof(REF_DBL));
	      RNS(xyz,"malloc xyz on master failed");
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
	  xyz=(REF_DBL*)malloc(3*n*sizeof(REF_DBL));
	  RNS(xyz,"malloc xyz on worker failed");
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

  chunk = MAX(10000, ntet/ref_mpi_n);

  ref_cell = ref_grid_tet(ref_grid);
  node_per = ref_cell_node_per(ref_cell);

  sent_c2n =(REF_INT *)malloc(node_per*chunk*sizeof(REF_INT));
  RNS(sent_c2n,"malloc failed");

  if ( ref_mpi_master )
    {
      long offset;
      offset = 4*7
	     + 8*3*nnode
             + 4*4*ntri
             + 4*5*nqua;
      fseek(file,offset,SEEK_SET);

      elements_to_send =(REF_INT *)malloc(ref_mpi_n*sizeof(REF_INT));
      RNS(elements_to_send,"malloc failed");

      c2n =(REF_INT *)malloc(node_per*chunk*sizeof(REF_INT));
      RNS(c2n,"malloc failed");

      ncell = 0;
      while ( ncell < ntet )
	{
	  section_size = MIN(chunk,ntet-ncell);
	  ncell += section_size;

	  for ( part = 0; part<ref_mpi_n ; part++ )
	    elements_to_send[part] = 0;
	  
	  for (cell=0;cell<section_size;cell++)
	    {
	      for (node=0;node<node_per;node++)
		{
		  RES(1, fread( &(c2n[node+node_per*cell]), 
				sizeof(REF_INT), 1, file ), "cn" );
		  SWAP_INT(c2n[node+node_per*cell]);
		  c2n[node+node_per*cell]--;
		  all_procs[node] = ref_part_implicit( nnode, ref_mpi_n, 
						       c2n[node+node_per*cell]);
		}
	      RSS( ref_sort_unique( node_per, all_procs, 
				    &cell_procs, unique_procs), "uniq" );
	      for (node=0;node<cell_procs;node++)
		elements_to_send[unique_procs[node]]++;
	    }
	  
	  part=0;
	  if ( elements_to_send[part] > 0 )
	    {
	      for (cell=0;cell<section_size;cell++)
		{
		  needcell = REF_FALSE;
		  for (node=0;node<node_per;node++)
		    needcell = needcell ||
		      ( part == ref_part_implicit( nnode, ref_mpi_n, 
						   c2n[node+node_per*cell]) );
		  if ( needcell )
		    {
		      RSS( ref_cell_add( ref_cell, &(c2n[node_per*cell]), 
					 &new_cell ), "add cell to off proc");
		    }
		}
	    } 

	  for ( part = 1; part<ref_mpi_n ; part++ )
	    if ( elements_to_send[part] > 0 ) 
	      {
		RSS( ref_mpi_send( &(elements_to_send[part]), 
				   1, REF_INT_TYPE, part ), "send" );
		send_size = 0;
		for (cell=0;cell<section_size;cell++)
		  {
		    needcell = REF_FALSE;
		    for (node=0;node<node_per;node++)
		      needcell = needcell ||
			( part == ref_part_implicit( nnode, ref_mpi_n, 
						     c2n[node+node_per*cell]) );
		    if ( needcell )
		      {
			for (node=0;node<node_per;node++)
			  sent_c2n[node+node_per*send_size] = 
			    c2n[node+node_per*cell];
			send_size++;
		      }
		
		  }
		RSS( ref_mpi_send( sent_c2n, node_per*send_size, 
				   REF_INT_TYPE, part ), "send" );
	      }
		
	}

      free(c2n);
      free(elements_to_send);

      /* signal we are done */
      for ( part = 1; part<ref_mpi_n ; part++ )
	RSS( ref_mpi_send( &end_of_message, 1, REF_INT_TYPE, part ), "send" );

    }  
  else
    {
      do {
	RSS( ref_mpi_recv( &elements_to_receive, 1, REF_INT_TYPE, 0 ), "recv" );
	printf("part %d, inbound %d\n",ref_mpi_id,elements_to_receive);
	if ( elements_to_receive > 0 )
	  {
	    RSS( ref_mpi_recv( sent_c2n, node_per*elements_to_receive, 
			       REF_INT_TYPE, 0 ), "send" );

	    for (cell=0;cell<elements_to_receive;cell++)
	      RSS( ref_cell_add( ref_cell, &(sent_c2n[node_per*cell]), 
				 &new_cell ), "add cell to off proc");
	  }
      } while ( elements_to_receive != end_of_message );
    }

  free(sent_c2n);

  if ( ref_mpi_master ) fclose(file);

  ref_node_n_global(ref_node) = nnode;

  each_ref_node_valid_node( ref_node, node )
    ref_node_part(ref_node,node) = 
    ref_part_implicit( nnode, ref_mpi_n, ref_node_global(ref_node,node));

  /* ghost xyz */

  RSS( ref_part_ghost_xyz( ref_grid ), "ghost xyz");

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

  a_size =(REF_INT *)malloc(ref_mpi_n*sizeof(REF_INT));
  RNS(a_size,"malloc failed");
  for ( part = 0; part<ref_mpi_n ; part++ )
    a_size[part] = 0;

  b_size =(REF_INT *)malloc(ref_mpi_n*sizeof(REF_INT));
  RNS(b_size,"malloc failed");
  for ( part = 0; part<ref_mpi_n ; part++ )
    b_size[part] = 0;

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      a_size[ref_node_part(ref_node,node)]++;

  RSS( ref_mpi_alltoall( a_size, b_size, REF_INT_TYPE ), "alltoall sizes");

  a_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    a_total += a_size[part];
  a_global =(REF_INT *)malloc(a_total*sizeof(REF_INT));
  RNS(a_global,"malloc failed");
  a_xyz =(REF_DBL *)malloc(a_total*3*sizeof(REF_DBL));
  RNS(a_xyz,"malloc failed");

  b_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    b_total += b_size[part];
  b_global =(REF_INT *)malloc(b_total*sizeof(REF_INT));
  RNS(b_global,"malloc failed");
  b_xyz =(REF_DBL *)malloc(b_total*3*sizeof(REF_DBL));
  RNS(b_xyz,"malloc failed");

  a_next =(REF_INT *)malloc(ref_mpi_n*sizeof(REF_INT));
  RNS(a_next,"malloc failed");
  a_next[part] = 0;

  a_next[0] = 0;
  for ( part = 1; part<ref_mpi_n ; part++ )
    a_next[part] = a_next[part-1]+a_size[part-1];

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      {
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
			  3, REF_INT_TYPE ), 
       "alltoallv global");

  for (node=0;node<a_total;node++)
    {
      RSS( ref_node_local( ref_node, a_global[node], &local ), "g2l");
      ref_node_xyz(ref_node,0,local) = a_xyz[0+3*node];
      ref_node_xyz(ref_node,1,local) = a_xyz[1+3*node];
      ref_node_xyz(ref_node,2,local) = a_xyz[2+3*node];
    }

  free(b_xyz);
  free(b_global);
  free(a_xyz);
  free(a_global);
  free(b_size);
  free(a_size);

  return REF_SUCCESS;  
}
