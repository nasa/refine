
#include <stdlib.h>
#include <stdio.h>

#include "ref_migrate.h"

REF_STATUS ref_migrate_create( REF_MIGRATE *ref_migrate_ptr, REF_GRID ref_grid )
{
  REF_MIGRATE ref_migrate;

  (*ref_migrate_ptr) = NULL;
  (*ref_migrate_ptr) = (REF_MIGRATE)malloc( sizeof(REF_MIGRATE_STRUCT) );
  RNS(*ref_migrate_ptr,"malloc ref_migrate NULL");

  ref_migrate = *ref_migrate_ptr;

  ref_migrate_grid(ref_migrate) = ref_grid;

  RSS( ref_edge_create( &(ref_migrate_edge( ref_migrate )), 
			ref_migrate_grid(ref_migrate) ), "create edge" );

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_free( REF_MIGRATE ref_migrate )
{
  if ( NULL == (void *)ref_migrate ) return REF_NULL;

  RSS( ref_edge_free( ref_migrate_edge( ref_migrate ) ), "free edge" );

  ref_cond_free( ref_migrate );

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_inspect( REF_MIGRATE ref_migrate )
{
  RSS(ref_grid_inspect( ref_migrate_grid(ref_migrate)) , "inspect grid");

  return REF_SUCCESS;
}

