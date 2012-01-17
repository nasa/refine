
#include <stdlib.h>
#include <stdio.h>

#include "ref_part.h"
#include "ref_mpi.h"
#include "ref_endian.h"

REF_STATUS ref_part_b8_ugrid( REF_GRID *ref_grid_ptr, char *filename )
{
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_INT nodes[REF_CELL_MAX_NODE_PER];
  REF_INT node, cell, new_cell;

  REF_GRID chunk;
  REF_CELL ref_cell;

  RSS( ref_grid_create( ref_grid_ptr ), "create grid");

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

  if ( ref_mpi_master )
    {
      long offset;
      offset = 4*7
	     + 8*3*nnode
             + 4*4*ntri
             + 4*5*nqua;
      fseek(file,offset,SEEK_SET);

      RSS( ref_grid_create( &chunk ), "create chunk");

      ref_cell = ref_grid_tet(chunk);
      for (cell=0;cell<ntri;cell++)
	{
	  for (node=0;node<ref_cell_node_per(ref_cell);node++)
	    {
	      RES(1, fread( &(nodes[node]), sizeof(REF_INT), 1, file ), "cn" );
	      SWAP_INT(nodes[node]);
	      nodes[node]--;
	    }
	  RSS( ref_cell_add(ref_cell, nodes, &new_cell ), "add cell");
	}
      RSS( ref_grid_free( chunk ), "free chunk");

   }  


  if ( ref_mpi_master ) fclose(file);

  return REF_SUCCESS;
}
