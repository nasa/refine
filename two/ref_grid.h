
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
  REF_NODE node;

  REF_CELL cell[4];

  REF_CELL tri;
  REF_CELL qua;

  REF_INT nedge;
};

REF_STATUS ref_grid_create( REF_GRID *ref_grid );
REF_STATUS ref_grid_free( REF_GRID ref_grid );

#define ref_grid_node(ref_grid) ((ref_grid)->node)

#define ref_grid_tet(ref_grid) ((ref_grid)->cell[0])
#define ref_grid_pyr(ref_grid) ((ref_grid)->cell[1])
#define ref_grid_pri(ref_grid) ((ref_grid)->cell[2])
#define ref_grid_hex(ref_grid) ((ref_grid)->cell[3])

#define ref_grid_tri(ref_grid) ((ref_grid)->tri)
#define ref_grid_qua(ref_grid) ((ref_grid)->qua)

#define ref_grid_nedge(ref_grid) ((ref_grid)->nedge)


#define each_ref_grid_ref_cell( ref_grid, group, ref_cell )		\
  for ( (group) = 0, (ref_cell) = (ref_grid)->cell[(group)] ;		\
	(group) < 4;							\
	(group)++  , (ref_cell) = (ref_grid)->cell[(group)] )

REF_STATUS ref_grid_inspect( REF_GRID ref_grid );

REF_STATUS ref_grid_make_edges( REF_GRID ref_grid );

END_C_DECLORATION

#endif /* REF_GRID_H */
