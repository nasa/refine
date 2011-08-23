
#ifndef REF_ADJ_H
#define REF_ADJ_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_ADJ_STRUCT REF_ADJ_STRUCT;
typedef REF_ADJ_STRUCT * REF_ADJ;
typedef struct REF_ADJ_ITEM_STRUCT REF_ADJ_ITEM_STRUCT;
typedef REF_ADJ_ITEM_STRUCT * REF_ADJ_ITEM;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_ADJ_STRUCT {
  REF_INT nnode, nitem;
  REF_INT *first;
  REF_ADJ_ITEM item;
  REF_INT blank;
};

struct REF_ADJ_ITEM_STRUCT {
  REF_INT next;
  REF_INT ref;
};

REF_STATUS ref_adj_create( REF_ADJ *ref_adj );
REF_STATUS ref_adj_free( REF_ADJ ref_adj );

#define ref_adj_nnode( ref_adj ) ((ref_adj)->nnode)
#define ref_adj_nitem( ref_adj ) ((ref_adj)->nitem)

#define ref_adj_first( ref_adj, node ) ( (ref_adj)->first[(node)] )
#define ref_adj_valid( item ) ( REF_EMPTY != (item) )
#define ref_adj_next( ref_adj, item ) ( (ref_adj)->item[(item)].next )
#define ref_adj_ref( ref_adj, item ) ( ref_adj_valid( item )?(ref_adj)->item[(item)].ref:REF_EMPTY )

#define ref_adj_for_loop( ref_adj, node, item, ref) \
  for ( (item) = ref_adj_first( ref_adj, node ),    \
	  (ref) = ref_adj_ref( ref_adj, item ) ;    \
	ref_adj_valid( item ) ;			    \
	(item) = ref_adj_next( ref_adj, item ),	    \
	  (ref) =  ref_adj_ref( ref_adj, item ) )

REF_STATUS ref_adj_inspect( REF_ADJ ref_adj );

REF_STATUS ref_adj_reg( REF_ADJ ref_adj, REF_INT node, REF_INT reference );


END_C_DECLORATION

#endif /* REF_ADJ_H */
