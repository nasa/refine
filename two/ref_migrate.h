
#ifndef REF_MIGRATE_H
#define REF_MIGRATE_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_MIGRATE_STRUCT REF_MIGRATE_STRUCT;
typedef REF_MIGRATE_STRUCT * REF_MIGRATE;
END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

struct REF_MIGRATE_STRUCT {
  REF_GRID grid;
  REF_INT *orig2;
  REF_INT n;
  REF_INT *global;
  REF_DBL *xyz;
  REF_DBL *weight;
};

#define ref_migrate_grid( ref_migrate ) ((ref_migrate)->grid)
#define ref_migrate_n( ref_migrate ) ((ref_migrate)->n)
#define ref_migrate_global( ref_migrate, node ) ((ref_migrate)->global[(node)])
#define ref_migrate_xyz( ref_migrate, ixyz, node ) \
  ((ref_migrate)->xyz[(ixyz)+3*(node)])
#define ref_migrate_weight( ref_migrate, node ) ((ref_migrate)->weight[(node)])

#define each_ref_migrate_node( ref_migrate, node )			\
  for ( (node) = 0 ;							\
	(node) < ref_migrate_n(ref_migrate);				\
	(node)++ )


REF_STATUS ref_migrate_create( REF_MIGRATE *ref_migrate, REF_GRID ref_grid );
REF_STATUS ref_migrate_free( REF_MIGRATE ref_migrate );

REF_STATUS ref_migrate_to_balance( REF_GRID ref_grid );

REF_STATUS ref_migrate_new_part( REF_GRID ref_grid );

REF_STATUS ref_migrate_shufflin( REF_GRID ref_grid );
REF_STATUS ref_migrate_shufflin_cell( REF_NODE ref_node, 
				      REF_CELL ref_cell );

END_C_DECLORATION

#endif /* REF_MIGRATE_H */
