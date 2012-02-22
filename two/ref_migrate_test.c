#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_metric.h"
#include "ref_sort.h"

#include "ref_migrate.h"

#include "ref_fixture.h"
#include "ref_export.h"
#include "ref_dict.h"

#include "ref_test.h"

static REF_STATUS set_up_tet_for_migrate( REF_MIGRATE *ref_migrate_ptr )
{
  REF_GRID ref_grid;

  TSS(ref_fixture_tet_grid( &ref_grid ), "tet");

  TSS(ref_migrate_create(ref_migrate_ptr,ref_grid),"create");

  return REF_SUCCESS;
}

static REF_STATUS tear_down( REF_MIGRATE ref_migrate )
{
  REF_GRID ref_grid;

  ref_grid = ref_migrate_grid(ref_migrate);

  TSS(ref_migrate_free(ref_migrate),"free");

  TSS( ref_grid_free(ref_grid),"free" );

  return REF_SUCCESS;
}

int main( void )
{

  { /* split tet in two, map 1 */
    REF_MIGRATE ref_migrate;
    REF_GRID ref_grid;

    TSS(set_up_tet_for_migrate(&ref_migrate),"set up");
    ref_grid = ref_migrate_grid(ref_migrate);

    TEIS(1, ref_cell_n(ref_grid_tet(ref_grid)),"two tet");
    TEIS(1, ref_cell_n(ref_grid_tri(ref_grid)),"two tri");

    TSS( tear_down( ref_migrate ), "tear down");
  }
  return 0;
}
