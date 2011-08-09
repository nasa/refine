
#ifndef CELL_H
#define CELL_H

#include "ref_def.h"

BEGIN_C_DECLORATION
typedef struct REF_CELL_STRUCT REF_CELL_STRUCT;
typedef REF_CELL_STRUCT * REF_CELL;
END_C_DECLORATION

BEGIN_C_DECLORATION
struct REF_CELL_STRUCT {
  REF_INT nodes, n, max;
  REF_INT *c2n
};

REF_STATUS ref_cell_create( REF_INT nodes, REF_CELL *ref_cell );
