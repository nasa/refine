
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_gather.h"

#include "ref_endian.h"
#include "ref_mpi.h"

REF_STATUS ref_gather_b8_ugrid( REF_GRID ref_grid, char *filename  )
{
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nnode;
  
  if ( ref_mpi_master )
    {
      file = fopen(filename,"w");
      if (NULL == (void *)file) printf("unable to open %s\n",filename);
      RNS(file, "unable to open file" );

      nnode = ref_node_n(ref_node);
      SWAP_INT(nnode);
      REIS(1, fwrite(&nnode,sizeof(REF_INT),1,file),"nnode");

    }

  if ( ref_mpi_master ) fclose(file);

  return REF_SUCCESS;
}

