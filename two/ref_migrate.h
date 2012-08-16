
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
};

#define ref_migrate_grid( ref_migrate ) ((ref_migrate)->grid)

REF_STATUS ref_migrate_create( REF_MIGRATE *ref_migrate, REF_GRID ref_grid );
REF_STATUS ref_migrate_free( REF_MIGRATE ref_migrate );

REF_STATUS ref_migrate_to_balance( REF_GRID ref_grid );

REF_STATUS ref_migrate_new_part( REF_GRID ref_grid );

REF_STATUS ref_migrate_shufflin( REF_GRID ref_grid );
REF_STATUS ref_migrate_shufflin_cell( REF_NODE ref_node, 
				      REF_CELL ref_cell );

END_C_DECLORATION

#endif /* REF_MIGRATE_H */
