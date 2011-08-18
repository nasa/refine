
#include <stdlib.h>
#include <stdio.h>

#include "ref_grid_import.h"

REF_STATUS ref_grid_export_vtk( REF_GRID ref_grid, char *filename  )
{
  FILE *file;
  REF_NODE ref_node;
  REF_INT node;
  REF_INT *o2n;

  ref_node = ref_grid->nodes;

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  fprintf(file,"# vtk DataFile Version 2.0\n");
  fprintf(file,"ref_grid_export_vtk\n");
  fprintf(file,"ASCII\n");

  RSS( ref_node_compact( ref_node, &o2n), "compact" );

  fprintf(file,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(file,"POINTS %d double\n",ref_node_n(ref_node));

  for ( node = 0; node < ref_node_max(ref_node); node++ )
    if ( REF_EMPTY != o2n[node] )
      fprintf(file, " %.16e %.16e %.16e\n",
	      ref_node_xyz(ref_node,0,node),
	      ref_node_xyz(ref_node,1,node),
	      ref_node_xyz(ref_node,2,node) ) ;

  free(o2n);

  fclose(file);

  return REF_SUCCESS;
}
