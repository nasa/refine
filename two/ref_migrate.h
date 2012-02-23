
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
  REF_EDGE edge;
};

#define ref_migrate_grid( ref_migrate ) ((ref_migrate)->grid)
#define ref_migrate_edge( ref_migrate ) ((ref_migrate)->edge)

REF_STATUS ref_migrate_create( REF_MIGRATE *ref_migrate, REF_GRID ref_grid );
REF_STATUS ref_migrate_free( REF_MIGRATE ref_migrate );

REF_STATUS ref_migrate_inspect( REF_MIGRATE ref_migrate );

int ref_migrate_number_of_vertices( void *void_ref_migrate, int *ierr );

END_C_DECLORATION

#endif /* REF_MIGRATE_H */
