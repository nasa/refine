
#ifndef REF_RECOVER_H
#define REF_RECOVER_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_RECOVER_STRUCT REF_RECOVER_STRUCT;
typedef REF_RECOVER_STRUCT * REF_RECOVER;
END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

struct REF_RECOVER_STRUCT {
  REF_GRID ref_grid;
  REF_INT n;
};

REF_STATUS ref_recover_create( REF_RECOVER *ref_recover, REF_GRID ref_grid );
REF_STATUS ref_recover_free( REF_RECOVER ref_recover );

#define ref_recover_n( ref_recover ) ((ref_recover)->n)
#define ref_recover_grid( ref_recover ) ((ref_recover)->ref_grid)

REF_STATUS ref_recover_enclosing_triangle( REF_RECOVER ref_recover, 
					   REF_INT node,
					   REF_INT *enclosing_cell,
					   REF_DBL *enclosing_bary);
REF_STATUS ref_recover_insert_twod( REF_RECOVER ref_recover, REF_DBL *xz,
				    REF_INT *node );
REF_STATUS ref_recover_opposite_node( REF_GRID ref_grid, 
				      REF_INT node0, REF_INT *node1 );
REF_STATUS ref_recover_edge_twod( REF_RECOVER ref_recover, 
				  REF_INT node0, REF_INT node1 );

REF_STATUS ref_recover_first_block_twod( REF_RECOVER ref_recover, 
					 REF_INT node0, REF_INT node1,
					 REF_INT *block1, REF_INT *block2 );

END_C_DECLORATION

#endif /* REF_RECOVER_H */
