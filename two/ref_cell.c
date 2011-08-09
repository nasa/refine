
#include <stdlib.h>
#include <stdio.h>
#include "ref_cell.h"

REF_STATUS ref_cell_create( REF_INT node_per, REF_CELL *ref_cell_ptr )
{
  *ref_cell_ptr = (REF_CELL)malloc( sizeof(REF_CELL_STRUCT) );
  RNS(*ref_cell_ptr,"malloc NULL");

  (*ref_cell_ptr)->node_per = node_per;
  (*ref_cell_ptr)->n = 0;
  (*ref_cell_ptr)->max = 0;
  (*ref_cell_ptr)->c2n = NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_free( REF_CELL ref_cell )
{
  ref_cond_free( ref_cell->c2n );
  return REF_SUCCESS;
}

REF_STATUS ref_cell_add( REF_CELL ref_cell, REF_INT *nodes )
{
  return REF_SUCCESS;
}
