
#ifndef REF_GEOM_H
#define REF_GEOM_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_GEOM_STRUCT REF_GEOM_STRUCT;
typedef REF_GEOM_STRUCT * REF_GEOM;
  END_C_DECLORATION

#include "ref_grid.h"
#include "ref_adj.h"

BEGIN_C_DECLORATION

#define REF_GEOM_NODE ( 0 )
#define REF_GEOM_EDGE ( 1 )
#define REF_GEOM_FACE ( 2 )

  struct REF_GEOM_STRUCT {
    REF_INT n, max;
    REF_INT blank;
    REF_INT *descr;
    REF_DBL *param;
    REF_ADJ ref_adj;
  };

#define ref_geom_n(ref_geom)     (( ref_geom )->n )
#define ref_geom_max(ref_geom)   (( ref_geom )->max )
#define ref_geom_blank(ref_geom) (( ref_geom )->blank )
#define ref_geom_adj(ref_geom)   (( ref_geom )->ref_adj )

#define ref_geom_descr(ref_geom,attribute,geom) \
  (( ref_geom )->descr[( attribute )+3*( geom )] )

#define ref_geom_type(ref_geom,geom) \
  (ref_geom_descr((ref_geom),0,(geom)))
#define ref_geom_id(ref_geom,geom) \
  (ref_geom_descr((ref_geom),1,(geom)))
#define ref_geom_node(ref_geom,geom) \
  (ref_geom_descr((ref_geom),2,(geom)))

#define ref_geom_param(ref_geom,dimension,geom) \
  (( ref_geom )->param[( dimension )+2*( geom )] )

#define each_ref_geom_edge( ref_geom, geom )			\
  for ( ( geom ) = 0;						\
        ( geom ) < ref_geom_max(ref_geom);			\
        ( geom )++ )						\
    if ( REF_GEOM_EDGE == ref_geom_type( ref_geom, geom ) )

#define each_ref_geom_face( ref_geom, geom )			\
  for ( ( geom ) = 0;						\
        ( geom ) < ref_geom_max(ref_geom);			\
        ( geom )++ )						\
    if ( REF_GEOM_FACE == ref_geom_type( ref_geom, geom ) )

REF_STATUS ref_geom_create( REF_GEOM *ref_geom );
REF_STATUS ref_geom_free( REF_GEOM ref_geom );

REF_STATUS ref_geom_deep_copy( REF_GEOM *ref_geom, REF_GEOM original );

REF_STATUS ref_geom_inspect( REF_GEOM ref_geom );

REF_STATUS ref_geom_add( REF_GEOM ref_geom, REF_INT node,
			 REF_INT type, REF_INT id,
			 REF_DBL *param );

REF_STATUS ref_geom_remove( REF_GEOM ref_geom, REF_INT node,
			    REF_INT type, REF_INT id);

REF_STATUS ref_geom_find( REF_GEOM ref_geom, REF_INT node,
			  REF_INT type, REF_INT id,
			  REF_INT *geom);

REF_STATUS ref_geom_tuv( REF_GEOM ref_geom, REF_INT node,
			 REF_INT type, REF_INT id,
			 REF_DBL *param );

REF_STATUS ref_geom_egads_fixture( char *filename );

REF_STATUS ref_geom_brep_from_egads( REF_GRID *ref_grid, char *filename );
REF_STATUS ref_geom_tetgen_volume( REF_GRID ref_grid );
REF_STATUS ref_geom_grid_from_egads( REF_GRID *ref_grid, char *filename );

REF_STATUS ref_geom_edge_tec_zone( REF_GRID ref_grid, REF_INT id, FILE *file );
REF_STATUS ref_geom_face_tec_zone( REF_GRID ref_grid, REF_INT id, FILE *file );
REF_STATUS ref_geom_tec( REF_GRID ref_grid, char *filename  );

END_C_DECLORATION

#endif /* REF_GEOM_H */
