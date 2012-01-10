
#ifndef REF_EXPORT_H
#define REF_EXPORT_H

#include "ref_grid.h"
#include "ref_node.h"
#include "ref_cell.h"
#include "ref_adj.h"
#include "ref_metric.h"

BEGIN_C_DECLORATION

REF_STATUS ref_export_vtk( REF_GRID ref_grid, char *filename );

REF_STATUS ref_export_tec( REF_GRID ref_grid, char *filename );
REF_STATUS ref_export_tec_surf( REF_GRID ref_grid, char *filename );

REF_STATUS ref_export_tec_surf_zone( REF_GRID ref_grid, FILE *file );
REF_STATUS ref_export_tec_vol_zone( REF_GRID ref_grid, FILE *file );

REF_STATUS ref_export_fgrid( REF_GRID ref_grid, char *filename );
REF_STATUS ref_export_ugrid( REF_GRID ref_grid, char *filename );
REF_STATUS ref_export_b8_ugrid( REF_GRID ref_grid, char *filename );

END_C_DECLORATION

#endif /* REF_EXPORT_H */
