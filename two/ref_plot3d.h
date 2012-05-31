
#ifndef REF_PLOT3D_H
#define REF_PLOT3D_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_PATCH_STRUCT REF_PATCH_STRUCT;
typedef REF_PATCH_STRUCT * REF_PATCH;

typedef struct REF_PLOT3D_STRUCT REF_PLOT3D_STRUCT;
typedef REF_PLOT3D_STRUCT * REF_PLOT3D;
END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION
struct REF_PLOT3D_STRUCT {
  REF_INT ngrid;
  REF_PATCH *patch;
};

struct REF_PATCH_STRUCT {
  REF_INT idim, jdim, kdim;
  REF_DBL *xyz;
};

REF_STATUS ref_plot3d_from_file( REF_PLOT3D *ref_plot3d, char *filename );
REF_STATUS ref_plot3d_free( REF_PLOT3D ref_plot3d );
REF_STATUS ref_patch_free( REF_PATCH ref_patch );

#define ref_plot3d_ngrid(ref_plot3d) ((ref_plot3d)->ngrid)

#define ref_patch_xyz(ref_patch,ixyz,i,j) \
  ((ref_patch)->xyz[ (ixyz) + 3 * ( (i) + (j)*(((ref_patch)->idim)) ) ] )

REF_STATUS ref_plot3d_tec( REF_PLOT3D ref_plot3d, char *filename );
REF_STATUS ref_plot3d_mate( REF_PLOT3D ref_plot3d, REF_GRID ref_grid );

REF_STATUS ref_patch_locate( REF_PATCH ref_patch, REF_DBL *xyz, REF_DBL *uv );
REF_STATUS ref_patch_xyz_at( REF_PATCH ref_patch, REF_DBL *uv, REF_DBL *xyz );
REF_STATUS ref_patch_dxyz_duv( REF_PATCH ref_patch, REF_DBL *uv, 
			       REF_DBL *dxyz_du, REF_DBL *dxyz_dv );

END_C_DECLORATION

#endif /* REF_PLOT3D_H */
