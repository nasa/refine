
#include <stdlib.h>
#include <stdio.h>

#include "ref_part.h"
#include "ref_mpi.h"
#include "ref_endian.h"

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
  REF_INT *elements_to_send;
  REF_INT node_per;
  REF_INT ncell, section_size;

  REF_INT chunk;

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
      for (node=0;node<ref_part_first( nnode, ref_mpi_n, 1 ); node++)
	{
	  RSS( ref_node_add(ref_node, node, &new_node ), "new_node");
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
	      ref_node_xyz( ref_node, 0, new_node ) = xyz[0+3*node];
	      ref_node_xyz( ref_node, 1, new_node ) = xyz[1+3*node];
	      ref_node_xyz( ref_node, 2, new_node ) = xyz[2+3*node];
	    }
	  free(xyz);
	}

    }


  if ( ref_mpi_master )
    {
      long offset;
      offset = 4*7
	     + 8*3*nnode
             + 4*4*ntri
             + 4*5*nqua;
      fseek(file,offset,SEEK_SET);

      elements_to_send =(REF_INT *)malloc(ref_mpi_n*sizeof(REF_INT));
      RNS(elements_to_send,"malloc elements_to_send on worker failed");

      ref_cell = ref_grid_tet(ref_grid);
      node_per = ref_cell_node_per(ref_cell);

      chunk = MAX(10000, ntet/ref_mpi_n);

      ncell = 0;
      while ( ncell < ntet )
	{
	  section_size = MIN(chunk,ntet-ncell);
	  ncell += section_size;
	  c2n =(REF_INT *)malloc(node_per*section_size*sizeof(REF_INT));
	  
	  for (cell=0;cell<section_size;cell++)
	    {
	      for (node=0;node<node_per;node++)
		{
		  RES(1, fread( &(c2n[node+node_per*cell]), 
				sizeof(REF_INT), 1, file ), "cn" );
		  SWAP_INT(c2n[node+node_per*cell]);
		  c2n[node+node_per*cell]--;
		}
	    }
	  free(c2n);
	}

      for ( part = 1; part<ref_mpi_n ; part++ )
	RSS( ref_mpi_send( &end_of_message, 1, REF_INT_TYPE, part ), "send" );

      free(elements_to_send);
    }  
  else
    {
      do {
	RSS( ref_mpi_recv( &elements_to_receive, 1, REF_INT_TYPE, 0 ), "recv" );
      } while ( elements_to_receive != end_of_message );
    }

  if ( ref_mpi_master ) fclose(file);

  return REF_SUCCESS;
}
