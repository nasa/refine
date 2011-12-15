
#ifndef REF_HEXDIV_H
#define REF_HEXDIV_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_HEXDIV_STRUCT REF_HEXDIV_STRUCT;
typedef REF_HEXDIV_STRUCT * REF_HEXDIV;
END_C_DECLORATION

#include "ref_grid.h"
#include "ref_face.h"

BEGIN_C_DECLORATION

struct REF_HEXDIV_STRUCT {
  REF_GRID grid;
  REF_FACE face;
  REF_INT *mark;
};

#define ref_hexdiv_grid( ref_hexdiv ) ((ref_hexdiv)->grid)
#define ref_hexdiv_face( ref_hexdiv ) ((ref_hexdiv)->face)

REF_STATUS ref_hexdiv_create( REF_HEXDIV *ref_hexdiv, REF_GRID ref_grid );
REF_STATUS ref_hexdiv_free( REF_HEXDIV ref_hexdiv );

#define ref_hexdiv_mark( ref_hexdiv, face )	\
  ((ref_hexdiv)->mark[face])

REF_STATUS ref_hexdiv_mark_to_split( REF_HEXDIV ref_hexdiv, 
				     REF_INT node0, REF_INT node1 );

REF_STATUS ref_hexdiv_mark_cell_edge_split( REF_HEXDIV ref_hexdiv, 
					    REF_INT cell, REF_INT cell_edge );

REF_STATUS ref_hexdiv_mark_relax( REF_HEXDIV ref_hexdiv );

END_C_DECLORATION

#endif /* REF_HEXDIV_H */
