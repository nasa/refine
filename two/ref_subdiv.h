
#ifndef REF_SUBDIV_H
#define REF_SUBDIV_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_SUBDIV_STRUCT REF_SUBDIV_STRUCT;
typedef REF_SUBDIV_STRUCT * REF_SUBDIV;
END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

struct REF_SUBDIV_STRUCT {
  REF_GRID ref_grid;
  REF_INT *mark;
};

#define ref_subdiv_grid(ref_subdiv) ((ref_subdiv)->ref_grid)

REF_STATUS ref_subdiv_create( REF_SUBDIV *ref_subdiv, REF_GRID ref_grid );
REF_STATUS ref_subdiv_free( REF_SUBDIV ref_subdiv );

END_C_DECLORATION

#endif /* REF_SUBDIV_H */
