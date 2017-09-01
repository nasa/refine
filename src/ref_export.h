
#ifndef REF_EXPORT_H
#define REF_EXPORT_H

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_export_by_extension( REF_GRID ref_grid, const char *filename );

REF_STATUS ref_export_vtk( REF_GRID ref_grid, const char *filename );

REF_STATUS ref_export_tec( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_tec_surf( REF_GRID ref_grid, const char *filename );

REF_STATUS ref_export_tec_surf_zone( REF_GRID ref_grid, FILE *file );
REF_STATUS ref_export_tec_vol_zone( REF_GRID ref_grid, FILE *file );

REF_STATUS ref_export_tec_int( REF_GRID ref_grid, REF_INT *scalar,
			       const char *filename );

REF_STATUS ref_export_tec_part( REF_GRID ref_grid, const char *root_filename );
REF_STATUS ref_export_metric_xyzdirlen( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_tec_metric_axis( REF_GRID ref_grid, const char *root_filename );
REF_STATUS ref_export_tec_metric_ellipse( REF_GRID ref_grid, 
					  const char *root_filename );
REF_STATUS ref_export_tec_metric_box( REF_GRID ref_grid, 
				      const char *root_filename,
				      REF_DBL *bounding_box );
REF_STATUS ref_export_tec_ratio( REF_GRID ref_grid, const char *root_filename );

REF_STATUS ref_export_smesh( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_fgrid( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_ugrid( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_lb8_ugrid( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_b8_ugrid( REF_GRID ref_grid, const char *filename );

REF_STATUS ref_export_cogsg( REF_GRID ref_grid, const char *filename );

REF_STATUS ref_export_c( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_eps( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_pdf( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_html( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_meshb( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_twod_meshb( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_twod_msh( REF_GRID ref_grid, const char *filename );

REF_STATUS ref_export_metric2d( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_twod_sol( REF_GRID ref_grid, const char *filename );

REF_STATUS ref_export_plt( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_export_plt_tet_zone( REF_GRID ref_grid, FILE *file );
REF_STATUS ref_export_plt_surf_zone( REF_GRID ref_grid, FILE *file );

REF_STATUS ref_export_faceid_range( REF_GRID ref_grid, 
				    REF_INT *min_faceid, REF_INT *max_faceid );
REF_STATUS ref_export_edgeid_range( REF_GRID ref_grid, 
				    REF_INT *min_edgeid, REF_INT *max_edgeid );

END_C_DECLORATION

#endif /* REF_EXPORT_H */
