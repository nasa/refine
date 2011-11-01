#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_adj.h"

#include "ref_subdiv.h"

#include "ref_test.h"

int main( void )
{
  REF_GRID ref_grid;
  REF_INT nodes[6] = {0,1,2,3,4,5};
  REF_INT cell;

  TSS(ref_grid_create(&ref_grid),"create");
  TSS(ref_cell_add(ref_grid_pri(ref_grid),nodes,&cell),"add prism");

  {
    REF_SUBDIV ref_subdiv;
    TSS(ref_subdiv_create(&ref_subdiv,ref_grid),"create");

    TES(0,ref_subdiv_mark(ref_subdiv,0),"init mark");
    TES(0,ref_subdiv_mark(ref_subdiv,6),"init mark");

    TSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0-1");

    TES(1,ref_subdiv_mark(ref_subdiv,0),"been marked");
    TES(0,ref_subdiv_mark(ref_subdiv,6),"not been marked");

    TSS(ref_subdiv_mark_relax(ref_subdiv),"relax");

    TES(1,ref_subdiv_mark(ref_subdiv,0),"been marked");
    TES(1,ref_subdiv_mark(ref_subdiv,6),"been relaxed");

    TSS(ref_subdiv_free(ref_subdiv),"free");
  }

  {
    REF_SUBDIV ref_subdiv;
    TSS(ref_subdiv_create(&ref_subdiv,ref_grid),"create");

    TSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0-1");
    TSS(ref_subdiv_mark_to_split(ref_subdiv,1,2),"mark edge 1-2");

    TES(0,ref_subdiv_mark(ref_subdiv,1),"no yet");
    TSS(ref_subdiv_mark_relax(ref_subdiv),"relax");
    TES(1,ref_subdiv_mark(ref_subdiv,1),"yet");

    TSS(ref_subdiv_free(ref_subdiv),"free");
  }

  {
    REF_SUBDIV ref_subdiv;
    TSS(ref_subdiv_create(&ref_subdiv,ref_grid),"create");
    TSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0-1");

    TSS(ref_subdiv_split(ref_subdiv),"split");
    TES(2, ref_cell_n(ref_grid_pri(ref_grid)),"two");
    TSS(ref_subdiv_free(ref_subdiv),"free");
  }

  TSS( ref_grid_free(ref_grid),"free" );

  return 0;
}
