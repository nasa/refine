
#ifndef REF_SUBDIV_H
#define REF_SUBDIV_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_SUBDIV_STRUCT REF_SUBDIV_STRUCT;
typedef REF_SUBDIV_STRUCT * REF_SUBDIV;
END_C_DECLORATION

#include "ref_grid.h"
#include "ref_adj.h"

BEGIN_C_DECLORATION

struct REF_SUBDIV_STRUCT {
  REF_GRID ref_grid;
  REF_ADJ adj;
  REF_INT *e2n;
  REF_INT *mark;
  REF_INT *node;
};

#define ref_subdiv_grid(ref_subdiv) ((ref_subdiv)->ref_grid)

REF_STATUS ref_subdiv_create( REF_SUBDIV *ref_subdiv, REF_GRID ref_grid );
REF_STATUS ref_subdiv_free( REF_SUBDIV ref_subdiv );

#define ref_subdiv_e2n( ref_subdiv, edge, node ) \
  ((ref_subdiv)->e2n[node+2*edge])

#define ref_subdiv_mark( ref_subdiv, edge )	\
  ((ref_subdiv)->mark[edge])

#define ref_subdiv_node( ref_subdiv, edge )	\
  ((ref_subdiv)->node[edge])

#define ref_subdiv_adj( ref_subdiv ) ((ref_subdiv)->adj)

REF_STATUS ref_subdiv_inspect( REF_SUBDIV ref_subdiv );

REF_STATUS ref_subdiv_edge_with( REF_SUBDIV ref_subdiv, 
				 REF_INT node0, REF_INT node1,
				 REF_INT *edge );

REF_STATUS ref_subdiv_mark_to_split( REF_SUBDIV ref_subdiv, 
				     REF_INT node0, REF_INT node1 );

REF_STATUS ref_subdiv_mark_relax( REF_SUBDIV ref_subdiv );

REF_STATUS ref_subdiv_new_node( REF_SUBDIV ref_subdiv );

REF_STATUS ref_subdiv_split( REF_SUBDIV ref_subdiv );

END_C_DECLORATION

#endif /* REF_SUBDIV_H */
