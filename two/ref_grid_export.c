
#include <stdlib.h>
#include <stdio.h>

#include "ref_grid_import.h"

#define VTK_TETRA      (10)
#define VTK_HEXAHEDRON (12)
#define VTK_WEDGE      (13)
#define VTK_PYRAMID    (14)


REF_STATUS ref_grid_export_vtk( REF_GRID ref_grid, char *filename  )
{
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n;
  REF_INT ncell,size;
  REF_INT *nodes;
  REF_INT node_per, cell;

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

  ncell = 0;
  ncell += ref_cell_n(ref_grid->cells[4]);
  ncell += ref_cell_n(ref_grid->cells[5]);
  ncell += ref_cell_n(ref_grid->cells[6]);
  ncell += ref_cell_n(ref_grid->cells[8]);
  size = 0;
  size += 5*ref_cell_n(ref_grid->cells[4]);
  size += 6*ref_cell_n(ref_grid->cells[5]);
  size += 7*ref_cell_n(ref_grid->cells[6]);
  size += 9*ref_cell_n(ref_grid->cells[8]);

  fprintf(file,"CELLS %d %d\n",ncell,size);

  node_per = 4;
  ref_cell = ref_grid->cells[node_per];
  nodes = (REF_INT *) malloc( node_per * sizeof(REF_INT) );
  for ( cell = 0 ; cell < ref_cell_max(ref_cell) ; cell++ )
    if ( ref_cell_valid( ref_cell, cell ) )
      {
	fprintf(file," %d",node_per);
	RSS(ref_cell_nodes( ref_cell, cell, nodes ), "cell nodes")
	for ( node = 0; node < node_per; node++ )
	  fprintf(file," %d",o2n[nodes[node]]+1);
	fprintf(file,"\n");
      }
  free(nodes);

  fprintf(file,"CELL_TYPES %d\n",ncell);
  node_per = 4;
  ref_cell = ref_grid->cells[node_per];
  for ( cell = 0 ; cell < ref_cell_n(ref_cell) ; cell++ )
    fprintf(file," %d\n",VTK_TETRA);


  free(o2n);

  fclose(file);

  return REF_SUCCESS;
}
