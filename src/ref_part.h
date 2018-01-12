
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#ifndef REF_PART_H
#define REF_PART_H

#include "ref_defs.h"

#include "ref_grid.h"
#include "ref_node.h"
#include "ref_mpi.h"

BEGIN_C_DECLORATION

#define ref_part_first( total_things, total_parts, part ) \
  (MIN(((((total_things)-1)/(total_parts))+1)*(part),(total_things)))

#define ref_part_implicit( total_things, total_parts, thing ) \
  ((thing) / ((((total_things)-1)/(total_parts))+1) )

REF_STATUS ref_part_by_extension( REF_GRID *ref_grid, 
				  REF_MPI ref_mpi, const char *filename );

REF_STATUS ref_part_b8_ugrid_cell( REF_CELL ref_cell, REF_INT ncell,
				   REF_NODE ref_node, REF_INT nnode,
				   FILE *file, 
				   long conn_offset,
				   long faceid_offset );
REF_STATUS ref_part_ghost_xyz( REF_GRID ref_grid );
REF_STATUS ref_part_ghost_int( REF_GRID ref_grid, REF_INT *scalar );

REF_STATUS ref_part_bamg_metric( REF_GRID ref_grid, const char *filename );
REF_STATUS ref_part_metric( REF_NODE ref_node, const char *filename );
REF_STATUS ref_part_scalar( REF_NODE ref_node, 
			    REF_DBL *scalar, const char *filename );

END_C_DECLORATION

#endif /* REF_PART_H */
