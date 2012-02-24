
#ifndef REF_MIGRATE_H
#define REF_MIGRATE_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_MIGRATE_STRUCT REF_MIGRATE_STRUCT;
typedef REF_MIGRATE_STRUCT * REF_MIGRATE;
END_C_DECLORATION

#include "ref_grid.h"
#include "ref_edge.h"

BEGIN_C_DECLORATION

struct REF_MIGRATE_STRUCT {
  REF_GRID grid;
  REF_INT *part;
};

#define ref_migrate_grid( ref_migrate ) ((ref_migrate)->grid)
#define ref_migrate_local( ref_migrate, node ) ((ref_migrate)->local[(node)])
#define ref_migrate_part( ref_migrate, node ) ((ref_migrate)->part[(node)])

REF_STATUS ref_migrate_create( REF_MIGRATE *ref_migrate, REF_GRID ref_grid );
REF_STATUS ref_migrate_free( REF_MIGRATE ref_migrate );

REF_STATUS ref_migrate_inspect( REF_MIGRATE ref_migrate );

REF_STATUS ref_migrate_part_viz( REF_MIGRATE ref_migrate );

END_C_DECLORATION

#endif /* REF_MIGRATE_H */
