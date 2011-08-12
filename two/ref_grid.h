
#ifndef REF_GRID_H
#define REF_GRID_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_GRID_STRUCT REF_GRID_STRUCT;
typedef REF_GRID_STRUCT * REF_GRID;
END_C_DECLORATION

#include "ref_node.h"
#include "ref_cell.h"

BEGIN_C_DECLORATION

struct REF_GRID_STRUCT {
  REF_NODE nodes;
  REF_CELL cells[9];
  REF_CELL faces[5];
};

REF_STATUS ref_grid_create( REF_GRID *ref_grid );
REF_STATUS ref_grid_free( REF_GRID ref_grid );

REF_STATUS ref_grid_import_ugrid( char *filename, REF_GRID *ref_grid );

END_C_DECLORATION

#endif /* REF_GRID_H */
