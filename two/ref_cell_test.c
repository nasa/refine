#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_cell.h"

int main( int argc, char *argv[] )
{
  
  REF_CELL ref_cell;

  RSS(ref_cell_create(4,&ref_cell),"create");
  RES(0,ref_cell_n(ref_cell),"should init zero cells");

  return 0;
}
