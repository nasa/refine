
#ifndef REF_GATHER_H
#define REF_GATHER_H

#include "ref_defs.h"

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"

#include "ref_mpi.h"

BEGIN_C_DECLORATION

#define ref_gather_blocking_frame( ref_grid, zone_title )				\
  RSS( ref_gather_tec_movie_frame( ref_grid, zone_title ), "movie frame" )

REF_STATUS ref_gather_plot( REF_GRID ref_grid, const char *filename );

REF_STATUS ref_gather_tec_movie_record_button( REF_BOOL on_or_off );
REF_STATUS ref_gather_tec_movie_frame( REF_GRID ref_grid,
				       const char *zone_title );

REF_STATUS ref_gather_tec_part( REF_GRID ref_grid, const char *filename );

REF_STATUS ref_gather_meshb( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_gather_b8_ugrid( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_gather_metric( REF_GRID ref_grid, const char *filename );

REF_STATUS ref_gather_ncell( REF_NODE ref_node, REF_CELL ref_cell, 
			     REF_INT *ncell );
REF_STATUS ref_gather_ngeom( REF_NODE ref_node, REF_GEOM ref_geom, 
			     REF_INT type, REF_INT *ngeom );

REF_STATUS ref_gather_node( REF_NODE ref_node,
			    REF_BOOL swap_endian, REF_BOOL has_id, FILE *file );
REF_STATUS ref_gather_node_tec_part( REF_NODE ref_node, FILE *file );
REF_STATUS ref_gather_node_metric( REF_NODE ref_node, FILE *file );

REF_STATUS ref_gather_geom( REF_NODE ref_node, REF_GEOM ref_geom, 
			    REF_INT type, FILE *file );
REF_STATUS ref_gather_cell( REF_NODE ref_node, REF_CELL ref_cell, 
			    REF_BOOL faceid_insted_of_c2n, REF_BOOL always_id,
			    REF_BOOL swap_endian,
			    REF_BOOL select_faceid, REF_INT faceid,
			    FILE *file );
REF_STATUS ref_gather_cell_tec( REF_NODE ref_node, REF_CELL ref_cell, 
				FILE *file );

END_C_DECLORATION

#endif /* REF_GATHER_H */
