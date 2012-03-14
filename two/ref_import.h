
#ifndef REF_IMPORT_H
#define REF_IMPORT_H

#include "ref_grid.h"
#include "ref_dict.h"

BEGIN_C_DECLORATION

REF_STATUS ref_import_examine_header( char *filename );

REF_STATUS ref_import_by_extension( REF_GRID *ref_grid, char *filename );

REF_STATUS ref_import_fgrid( REF_GRID *ref_grid, char *filename );
REF_STATUS ref_import_ugrid( REF_GRID *ref_grid, char *filename );
REF_STATUS ref_import_b8_ugrid( REF_GRID *ref_grid, char *filename );
REF_STATUS ref_import_r8_ugrid( REF_GRID *ref_grid, char *filename );
REF_STATUS ref_import_msh( REF_GRID *ref_grid, char *filename );

REF_STATUS ref_import_mapbc( REF_DICT *ref_dict, char *filename );

END_C_DECLORATION

#endif /* REF_IMPORT_H */
