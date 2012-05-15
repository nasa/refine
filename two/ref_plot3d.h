
#ifndef REF_PLOT3D_H
#define REF_PLOT3D_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_PLOT3D_STRUCT REF_PLOT3D_STRUCT;
typedef REF_PLOT3D_STRUCT * REF_PLOT3D;

struct REF_PLOT3D_STRUCT {
  REF_INT ngrid;
  REF_INT *idim, *jdim, *kdim;
};

REF_STATUS ref_plot3d_from_file( REF_PLOT3D *ref_plot3d, char *filename );

#define ref_plot3d_ngrid(ref_plot3d) ((ref_plot3d)->ngrid)

END_C_DECLORATION

#endif /* REF_PLOT3D_H */
