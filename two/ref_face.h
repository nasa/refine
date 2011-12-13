
#ifndef REF_FACE_H
#define REF_FACE_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_FACE_STRUCT REF_FACE_STRUCT;
typedef REF_FACE_STRUCT * REF_FACE;
END_C_DECLORATION

#include "ref_grid.h"
#include "ref_adj.h"

BEGIN_C_DECLORATION

struct REF_FACE_STRUCT {
  REF_INT n, max;
  REF_INT *f2n;
  REF_ADJ adj;
};

REF_STATUS ref_face_create( REF_FACE *ref_face, REF_GRID ref_grid );
REF_STATUS ref_face_free( REF_FACE ref_face );

#define ref_face_n(ref_face) ((ref_face)->n)
#define ref_face_max(ref_face) ((ref_face)->max)

#define ref_face_f2n(ref_face,node,face) ((ref_face)->f2n[(node)+4*(face)])

#define ref_face_adj(ref_face) ((ref_face)->adj)

REF_STATUS ref_face_inspect( REF_FACE ref_face );

REF_STATUS ref_face_with( REF_FACE ref_face, REF_INT *nodes, REF_INT *face );
REF_STATUS ref_face_spanning( REF_FACE ref_face, REF_INT node0, REF_INT node1, 
			      REF_INT *face );

REF_STATUS ref_face_add_uniquely( REF_FACE ref_face, REF_INT *nodes );

REF_STATUS ref_face_normal( REF_DBL *xyz0, REF_DBL *xyz1, 
			    REF_DBL *xyz2, REF_DBL *xyz3, 
			    REF_DBL *normal );
END_C_DECLORATION

#endif /* REF_FACE_H */
