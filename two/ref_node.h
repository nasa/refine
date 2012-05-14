
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
  REF_DBL *real;
  REF_INT naux;
  REF_DBL *aux;
  REF_LIST unused_global_list;
  REF_INT old_n_global, new_n_global;
};

#define REF_NODE_REAL_PER (9) /* x,y,z, m[6] */

#define ref_node_n(ref_node) ((ref_node)->n)
#define ref_node_max(ref_node) ((ref_node)->max)

#define ref_node_n_global(ref_node) ((ref_node)->old_n_global)

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

#define ref_node_xyz(ref_node,ixyz,node) \
  ((ref_node)->real[(ixyz)+REF_NODE_REAL_PER*(node)])
#define ref_node_xyz_ptr(ref_node,node)			\
  (&((ref_node)->real[REF_NODE_REAL_PER*(node)]))

#define ref_node_metric(ref_node,im,node)		\
  ((ref_node)->real[(im+3)+REF_NODE_REAL_PER*(node)])
#define ref_node_metric_ptr(ref_node,node)		\
  (&((ref_node)->real[3+REF_NODE_REAL_PER*(node)]))

#define ref_node_real(ref_node,ireal,node)		\
  ((ref_node)->real[(ireal)+REF_NODE_REAL_PER*(node)])

#define ref_node_part(ref_node,node) ((ref_node)->part[(node)])

#define ref_node_naux(ref_node) ((ref_node)->naux)
#define ref_node_aux(ref_node,iaux,node)		\
  ((ref_node)->aux[(iaux)+ref_node_naux(ref_node)*(node)])

REF_STATUS ref_node_create( REF_NODE *ref_node );
REF_STATUS ref_node_free( REF_NODE ref_node );
REF_STATUS ref_node_inspect( REF_NODE ref_node );
REF_STATUS ref_node_location( REF_NODE ref_node, REF_INT node );
REF_STATUS ref_node_tattle_global( REF_NODE ref_node, REF_INT global );

REF_STATUS ref_node_local( REF_NODE ref_node, REF_INT global, REF_INT *node );

REF_STATUS ref_node_initialize_n_global(  REF_NODE ref_node, REF_INT n_global );
REF_STATUS ref_node_next_global( REF_NODE ref_node, REF_INT *global );

REF_STATUS ref_node_synchronize_globals( REF_NODE ref_node );
REF_STATUS ref_node_shift_new_globals( REF_NODE ref_node );
REF_STATUS ref_node_eliminate_unused_globals( REF_NODE ref_node );

REF_STATUS ref_node_add( REF_NODE ref_node, REF_INT global, REF_INT *node );
REF_STATUS ref_node_add_many( REF_NODE ref_node, REF_INT n, REF_INT *global );

REF_STATUS ref_node_remove( REF_NODE ref_node, REF_INT node );
REF_STATUS ref_node_remove_without_global( REF_NODE ref_node, REF_INT node );
REF_STATUS ref_node_rebuild_sorted_global( REF_NODE ref_node );

REF_STATUS ref_node_compact( REF_NODE ref_node, REF_INT **o2n, REF_INT **n2o );

REF_STATUS ref_node_ghost_real( REF_NODE ref_node );
REF_STATUS ref_node_ghost_int( REF_NODE ref_node, REF_INT *scalar );

REF_STATUS ref_node_ratio( REF_NODE ref_node, REF_INT node0, REF_INT node1, 
			   REF_DBL *ratio );
REF_STATUS ref_node_tri_normal( REF_NODE ref_node, 
				REF_INT *nodes,  
				REF_DBL *normal );
REF_STATUS ref_node_tet_vol( REF_NODE ref_node, 
			     REF_INT *nodes,  
			     REF_DBL *volume );
REF_STATUS ref_node_tet_quality( REF_NODE ref_node, 
				 REF_INT *nodes,  
				 REF_DBL *quality );

REF_STATUS ref_node_interpolate_edge( REF_NODE ref_node, 
				      REF_INT node0, REF_INT node1, 
				      REF_INT new_node );
REF_STATUS ref_node_resize_aux( REF_NODE ref_node );

END_C_DECLORATION

#endif /* REF_NODE_H */
