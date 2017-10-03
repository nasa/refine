
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

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

REF_STATUS ref_adj_deep_copy( REF_ADJ *ref_adj, REF_ADJ original );

#define ref_adj_nnode( ref_adj ) (( ref_adj )->nnode )
#define ref_adj_nitem( ref_adj ) (( ref_adj )->nitem )
#define ref_adj_blank( ref_adj ) (( ref_adj )->blank )

#define ref_adj_first( ref_adj, node )               \
  ( ( node )>=0&&( node )<ref_adj_nnode( ref_adj ) ? \
    ( ref_adj )->first[( node )] : REF_EMPTY )
#define ref_adj_valid( item ) ( REF_EMPTY != ( item ) )
#define ref_adj_item_next( ref_adj, item_arg ) \
  ( ( ref_adj )->item[( item_arg )].next )
#define ref_adj_item_ref( ref_adj, item_arg ) \
  ( ( ref_adj )->item[( item_arg )].ref )

#define ref_adj_safe_ref( ref_adj, item_arg ) \
  ( ref_adj_valid( item_arg ) ? ( ref_adj )->item[( item_arg )].ref : REF_EMPTY )

#define ref_adj_empty( ref_adj, node ) \
  ( !ref_adj_valid( ref_adj_first( ref_adj, node ) ) )

#define each_ref_adj_node_item_with_ref( ref_adj, node, item, ref) \
  for ( ( item ) = ref_adj_first( ref_adj, node ),                 \
        ( ref ) = ref_adj_safe_ref( ref_adj, item );               \
        ref_adj_valid( item );                                     \
        ( item ) = ref_adj_item_next( ref_adj, item ),             \
        ( ref ) =  ref_adj_safe_ref( ref_adj, item ) )

REF_STATUS ref_adj_inspect( REF_ADJ ref_adj );
REF_STATUS ref_adj_node_inspect( REF_ADJ ref_adj, REF_INT node );

REF_STATUS ref_adj_add( REF_ADJ ref_adj, REF_INT node, REF_INT reference );
REF_STATUS ref_adj_remove( REF_ADJ ref_adj, REF_INT node, REF_INT reference );

REF_STATUS ref_adj_add_uniquely( REF_ADJ ref_adj,
                                 REF_INT node, REF_INT reference );
REF_STATUS ref_adj_degree( REF_ADJ ref_adj,
                           REF_INT node, REF_INT *degree );

END_C_DECLORATION

#endif /* REF_ADJ_H */
