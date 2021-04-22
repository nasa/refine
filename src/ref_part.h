
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
#include "ref_mpi.h"
#include "ref_node.h"

BEGIN_C_DECLORATION

/* find the size of large parts: integer divide, round up */
#define ref_part_large_part_size(total_things, total_parts) \
  (((total_things) + (total_parts)-1) / (total_parts))

/* find the size of small parts: integer divide, round off */
#define ref_part_small_part_size(total_things, total_parts) \
  (((total_things)-1) / (total_parts))

/* find the number of large parts: remainder after filling small parts */
#define ref_part_n_large_part(total_things, total_parts) \
  ((total_things) -                                      \
   (total_parts)*ref_part_small_part_size(total_things, total_parts))

/* first thing index on a part, valid for 0 to nparts (returns total_things) */
#define ref_part_first(total_things, total_parts, part)                \
  ((part) < ref_part_n_large_part(total_things, total_parts)           \
       ? (part)*ref_part_large_part_size(total_things, total_parts)    \
       : ((part)-ref_part_n_large_part(total_things, total_parts)) *   \
                 ref_part_small_part_size(total_things, total_parts) + \
             (ref_part_n_large_part(total_things, total_parts) *       \
              ref_part_large_part_size(total_things, total_parts)))

#define ref_part_large_implicit(total_things, total_parts, thing) \
  ((thing) / ref_part_large_part_size(total_things, total_parts))

#define ref_part_total_large(total_things, total_parts) \
  (ref_part_n_large_part(total_things, total_parts) *   \
   ref_part_large_part_size(total_things, total_parts))

#define ref_part_small_implicit(total_things, total_parts, thing) \
  (((thing)-ref_part_total_large(total_things, total_parts)) /    \
       ref_part_small_part_size(total_things, total_parts) +      \
   ref_part_n_large_part(total_things, total_parts))

/* part id for a thing */
#define ref_part_implicit(total_things, total_parts, thing)                 \
  (REF_INT)((thing) / ref_part_large_part_size(total_things, total_parts) < \
                    ref_part_n_large_part(total_things, total_parts)        \
                ? ref_part_large_implicit(total_things, total_parts, thing) \
                : ref_part_small_implicit(total_things, total_parts, thing))

REF_STATUS ref_part_by_extension(REF_GRID *ref_grid, REF_MPI ref_mpi,
                                 const char *filename);

REF_STATUS ref_part_cad_data(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_part_cad_association(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_part_cad_discrete_edge(REF_GRID ref_grid, const char *filename);

REF_STATUS ref_part_bamg_metric(REF_GRID ref_grid, const char *filename);
REF_STATUS ref_part_metric(REF_NODE ref_node, const char *filename);
REF_STATUS ref_part_scalar(REF_GRID ref_grid, REF_INT *ldim, REF_DBL **scalar,
                           const char *filename);

END_C_DECLORATION

#endif /* REF_PART_H */
