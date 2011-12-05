#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid_export.h"
#include "ref_test.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node, cell;
  REF_INT nodes[8];

  if (argc>1) {printf("%s ignored\n",argv[0]);}

  TSS(ref_grid_create( &ref_grid ), "create" );
  ref_node = ref_grid_node(ref_grid);

  TSS(ref_node_add(ref_node,0,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  TSS(ref_node_add(ref_node,1,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  TSS(ref_node_add(ref_node,2,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  TSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  nodes[0]=0;
  nodes[1]=1;
  nodes[2]=2;
  nodes[3]=3;
  ref_cell = ref_grid_tet(ref_grid);
  TSS( ref_cell_add( ref_cell, nodes, &cell ), "add cell");

  TSS(ref_grid_export_vtk( ref_grid, "ref_grid_export_test.vtk" ),"export" );
  TSS(ref_grid_free(ref_grid),"free");

  return 0;
}
