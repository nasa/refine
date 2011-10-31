
#include <stdlib.h>
#include <stdio.h>

#include "ref_subdiv.h"

REF_STATUS ref_subdiv_create( REF_SUBDIV *ref_subdiv_ptr, REF_GRID ref_grid )
{
  REF_SUBDIV ref_subdiv;
  (*ref_subdiv_ptr) = NULL;
  (*ref_subdiv_ptr) = (REF_SUBDIV)malloc( sizeof(REF_SUBDIV_STRUCT) );
  RNS(*ref_subdiv_ptr,"malloc ref_subdiv NULL");

  ref_subdiv = *ref_subdiv_ptr;

  ref_subdiv_grid(ref_subdiv) = ref_grid;

  RSS( ref_grid_make_edges( ref_grid ), "edges");

  ref_subdiv->mark = (REF_INT *)malloc( ref_grid_nedge(ref_grid) 
					* sizeof(REF_INT));
  RNS(ref_subdiv->mark,"malloc global NULL");

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_free( REF_SUBDIV ref_subdiv )
{
  if ( NULL == (void *)ref_subdiv ) return REF_NULL;

  free( ref_subdiv->mark );

  ref_cond_free( ref_subdiv );

  return REF_SUCCESS;
}
