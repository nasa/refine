
#ifndef REF_NODE_H
#define REF_NODE_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_NODE_STRUCT REF_NODE_STRUCT;
typedef REF_NODE_STRUCT * REF_NODE;
END_C_DECLORATION

#include "ref_list.h"

BEGIN_C_DECLORATION

struct REF_NODE_STRUCT {
  REF_INT n, max;
  REF_INT blank;
  REF_INT *global;
  REF_INT *sorted_global;
  REF_INT *sorted_local;
  REF_INT *part;
  REF_DBL *xyz;
  REF_LIST unused_global_list;
  REF_INT old_n_global, new_n_global;
};

#define ref_node_n(ref_node) ((ref_node)->n)
#define ref_node_max(ref_node) ((ref_node)->max)

#define ref_node_valid(ref_node,node) \
  ( (node) > -1 && (node) < ref_node_max(ref_node) && \
    (ref_node)->global[(node)] >= 0 )

#define ref_node_global(ref_node,node)			\
  ( ( (node) > -1 && (node) < ref_node_max(ref_node) && \
      (ref_node)->global[(node)] >= 0) ?		\
    (ref_node)->global[(node)]:REF_EMPTY )

#define each_ref_node_valid_node( ref_node, node )			\
  for ( (node) = 0 ;							\
	(node) < ref_node_max(ref_node);				\
	(node)++ )							\
    if ( ref_node_valid( ref_node, node ) )

#define ref_node_xyz(ref_node,ixyz,node) ((ref_node)->xyz[(ixyz)+3*(node)])

#define ref_node_part(ref_node,node) ((ref_node)->part[(node)])

REF_STATUS ref_node_create( REF_NODE *ref_node );
REF_STATUS ref_node_free( REF_NODE ref_node );
REF_STATUS ref_node_inspect( REF_NODE ref_node );
REF_STATUS ref_node_location( REF_NODE ref_node, REF_INT node );

REF_STATUS ref_node_local( REF_NODE ref_node, REF_INT global, REF_INT *node );

REF_STATUS ref_node_initialize_n_global(  REF_NODE ref_node, REF_INT n_global );
REF_STATUS ref_node_next_global( REF_NODE ref_node, REF_INT *global );
REF_STATUS ref_node_sync_new_globals( REF_NODE ref_node );

REF_STATUS ref_node_add( REF_NODE ref_node, REF_INT global, REF_INT *node );
REF_STATUS ref_node_add_many( REF_NODE ref_node, REF_INT n, REF_INT *global );

REF_STATUS ref_node_remove( REF_NODE ref_node, REF_INT node );
REF_STATUS ref_node_remove_without_global( REF_NODE ref_node, REF_INT node );
REF_STATUS ref_node_rebuild_sorted_global( REF_NODE ref_node );

REF_STATUS ref_node_compact( REF_NODE ref_node, REF_INT *o2n[] );

END_C_DECLORATION

#endif /* REF_NODE_H */
