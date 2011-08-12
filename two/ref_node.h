
#ifndef REF_NODE_H
#define REF_NODE_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_NODE_STRUCT REF_NODE_STRUCT;
typedef REF_NODE_STRUCT * REF_NODE;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_NODE_STRUCT {
  REF_INT n, max;
  REF_INT *global;
};

REF_STATUS ref_node_create( REF_NODE *ref_node );
REF_STATUS ref_node_free( REF_NODE ref_node );
REF_STATUS ref_node_inspect( REF_NODE ref_node );

#define ref_node_n(ref_node) ((ref_node)->n)
#define ref_node_max(ref_node) ((ref_node)->max)

#define ref_node_valid(ref_node,node) \
  ( (node) > -1 && (node) < ref_node_max(ref_node) )

#define ref_node_global(ref_node,node) \
  ( ( (node) > -1 && (node) < ref_node_max(ref_node) ) ? \
           (ref_node)->global[(node)]:REF_EMPTY )

REF_STATUS ref_node_add( REF_NODE ref_node, REF_INT global, REF_INT *node );

END_C_DECLORATION

#endif /* REF_NODE_H */
