
#ifndef REF_CELL_H
#define REF_CELL_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_CELL_STRUCT REF_CELL_STRUCT;
typedef REF_CELL_STRUCT * REF_CELL;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_CELL_STRUCT {
  REF_INT nodes, n, max;
  REF_INT *c2n;
};

REF_STATUS ref_cell_create( REF_INT nodes, REF_CELL *ref_cell );

#define ref_cell_n(ref_cell) ((ref_cell)->n)

END_C_DECLORATION

#endif /* REF_CELL_H */
