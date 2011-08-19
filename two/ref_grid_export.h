
#ifndef REF_GRID_EXPORT_H
#define REF_GRID_EXPORT_H

#include "ref_grid.h"
#include "ref_node.h"
#include "ref_cell.h"

BEGIN_C_DECLORATION

REF_STATUS ref_grid_export_vtk( REF_GRID ref_grid, char *filename );
REF_STATUS ref_grid_export_fgrid( REF_GRID ref_grid, char *filename );

END_C_DECLORATION

#endif /* REF_GRID_EXPORT_H */
