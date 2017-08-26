
#ifndef REF_CAVITY_H
#define REF_CAVITY_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_CAVITY_STRUCT REF_CAVITY_STRUCT;
typedef REF_CAVITY_STRUCT * REF_CAVITY;
END_C_DECLORATION

#include "ref_grid.h"
#include "ref_list.h"

BEGIN_C_DECLORATION

struct REF_CAVITY_STRUCT {
  REF_INT node_per;
  REF_INT n;
  REF_INT max;
  REF_INT blank;
  REF_INT *f2n;
  REF_LIST ref_list;
};

REF_STATUS ref_cavity_create( REF_CAVITY *ref_cavity, REF_INT node_per );
REF_STATUS ref_cavity_free( REF_CAVITY ref_cavity );
REF_STATUS ref_cavity_inspect( REF_CAVITY ref_cavity );

#define ref_cavity_n( ref_cavity ) (( ref_cavity )->n )
#define ref_cavity_node_per( ref_cavity ) (( ref_cavity )->node_per )

#define ref_cavity_f2n(ref_cavity,node,cavity) \
  (( ref_cavity )->f2n[( node )+ref_cavity_node_per(ref_cavity)*( cavity )] )

#define ref_cavity_max( ref_cavity ) (( ref_cavity )->max )
#define ref_cavity_blank( ref_cavity ) (( ref_cavity )->blank )

#define ref_cavity_list( ref_cavity ) (( ref_cavity )->ref_list )

#define ref_cavity_valid(ref_cavity,face)                    \
  ( ( face ) >=0 && ( face ) < ref_cavity_max(ref_cavity) && \
    REF_EMPTY != ref_cavity_f2n(ref_cavity,0,face) )

#define each_ref_cavity_valid_face( ref_cavity, face ) \
  for ( ( face ) = 0;                                  \
        ( face ) < ref_cavity_max(ref_cavity);         \
        ( face )++ )                                   \
    if ( ref_cavity_valid( ref_cavity, face ) )

#define each_ref_cavity_face_node( ref_cavity, face_node ) \
  for ( ( face_node ) = 0;                                 \
        ( face_node ) < ref_cavity_node_per(ref_cavity);   \
        ( face_node )++ )

REF_STATUS ref_cavity_insert( REF_CAVITY ref_cavity, REF_INT *nodes );
REF_STATUS ref_cavity_find( REF_CAVITY ref_cavity, REF_INT *nodes,
                            REF_INT *found_face, REF_BOOL *reversed);

REF_STATUS ref_cavity_add_tet( REF_CAVITY ref_cavity,
                               REF_GRID ref_grid, REF_INT tet );
REF_STATUS ref_cavity_rm_tet( REF_CAVITY ref_cavity,
                              REF_GRID ref_grid, REF_INT tet );
REF_STATUS ref_cavity_replace_tet( REF_CAVITY ref_cavity,
                                   REF_GRID ref_grid, REF_INT node );

REF_STATUS ref_cavity_add_tri( REF_CAVITY ref_cavity,
                               REF_GRID ref_grid, REF_INT tri );
REF_STATUS ref_cavity_rm_tri( REF_CAVITY ref_cavity,
                              REF_GRID ref_grid, REF_INT tri );
REF_STATUS ref_cavity_add_ball( REF_CAVITY ref_cavity,
                                REF_GRID ref_grid, REF_INT node );

REF_STATUS ref_cavity_replace_tri( REF_CAVITY ref_cavity,
                                   REF_GRID ref_grid,
                                   REF_INT node, REF_INT clone );

REF_STATUS ref_cavity_visible( REF_CAVITY ref_cavity,
                               REF_NODE ref_node, REF_INT node, REF_INT face,
                               REF_BOOL *visible );
REF_STATUS ref_cavity_enlarge_visible( REF_CAVITY ref_cavity,
                                       REF_GRID ref_grid, REF_INT node );
REF_STATUS ref_cavity_shrink_visible( REF_CAVITY ref_cavity,
                                      REF_GRID ref_grid, REF_INT node );
REF_STATUS ref_cavity_make_visible( REF_CAVITY ref_cavity,
                                    REF_GRID ref_grid, REF_INT node );

REF_STATUS ref_cavity_enlarge_metric( REF_CAVITY ref_cavity,
                                      REF_GRID ref_grid, REF_INT node );

REF_STATUS ref_cavity_enlarge_face( REF_CAVITY ref_cavity,
                                    REF_GRID ref_grid, REF_INT face );
REF_STATUS ref_cavity_shrink_face( REF_CAVITY ref_cavity,
                                   REF_GRID ref_grid, REF_INT face );

REF_STATUS ref_cavity_tet_quality( REF_GRID ref_grid );
REF_STATUS ref_cavity_twod_pass( REF_GRID ref_grid );
REF_STATUS ref_cavity_tec( REF_CAVITY ref_cavity, REF_GRID ref_grid,
                           REF_INT node, const char *filename );

REF_STATUS ref_cavity_change( REF_CAVITY ref_cavity, REF_GRID ref_grid,
                              REF_INT node );



END_C_DECLORATION

#endif /* REF_CAVITY_H */
