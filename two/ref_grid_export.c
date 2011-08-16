
#include <stdlib.h>
#include <stdio.h>

#include "ref_grid_import.h"

REF_STATUS ref_grid_export_vtk( REF_GRID ref_grid, char *filename  )
{
  FILE *file;

  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  fprintf(file,"# vtk DataFile Version 2.0\n");
  fprintf(file,"ref_grid_export_vtk\n");
  fprintf(file,"ASCII\n");

  fprintf(file,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(file,"POINTS %d\n",ref_node_n(ref_grid->nodes));

  fclose(file);

  return REF_SUCCESS;
}
