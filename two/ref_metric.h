
#ifndef REF_METRIC_H
#define REF_METRIC_H

#include "ref_defs.h"

#include "ref_grid.h"
#include "ref_node.h"

BEGIN_C_DECLORATION

REF_STATUS ref_metric_show( REF_DBL *metric );
REF_STATUS ref_metric_inspect( REF_NODE ref_node );

REF_STATUS ref_metric_from_node( REF_DBL *metric, REF_NODE ref_node );
REF_STATUS ref_metric_to_node( REF_DBL *metric, REF_NODE ref_node );

REF_STATUS ref_metric_unit_node( REF_NODE ref_node );
REF_STATUS ref_metric_olympic_node( REF_NODE ref_node, REF_DBL h );
REF_STATUS ref_metric_polar2d_node( REF_NODE ref_node );
REF_STATUS ref_metric_ugawg_node( REF_NODE ref_node, REF_INT version );
REF_STATUS ref_metric_masabl_node( REF_NODE ref_node );
REF_STATUS ref_metric_twod_node( REF_NODE ref_node );

REF_STATUS ref_metric_interpolate( REF_GRID ref_grid, REF_GRID parent  );

REF_STATUS ref_metric_gradation( REF_GRID ref_grid, REF_DBL r );

REF_STATUS ref_metric_sanitize( REF_GRID ref_grid );
REF_STATUS ref_metric_sanitize_threed( REF_GRID ref_grid );
REF_STATUS ref_metric_sanitize_twod( REF_GRID ref_grid );

REF_STATUS ref_metric_from_curvature( REF_DBL *metric, REF_GRID ref_grid );

REF_STATUS ref_metric_imply_from( REF_DBL *metric, REF_GRID ref_grid );
REF_STATUS ref_metric_imply_non_tet( REF_DBL *metric, REF_GRID ref_grid );

REF_STATUS ref_metric_smr( REF_DBL *metric0, REF_DBL *metric1, REF_DBL *metric, 
			   REF_GRID ref_grid );

END_C_DECLORATION

#endif /* REF_METRIC_H */
