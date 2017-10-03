
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_import.h"
#include "ref_export.h"
#include "ref_validation.h"

#include "ref_adj.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_matrix.h"

#include "ref_cell.h"

#include "ref_face.h"
#include "ref_sort.h"
#include "ref_dict.h"
#include "ref_mpi.h"
#include "ref_edge.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc>1) 
    {
      printf("validating\n");

      printf("reading %s\n",argv[1]);
      RSS(ref_import_by_extension( &ref_grid, argv[1] ),"from ugrid");
      printf("complete.\n");
      
      RSS(ref_grid_inspect( ref_grid ), "inspection");

      printf("validate.\n");
      RSS( ref_validation_all( ref_grid ), "invalid grid" );

      printf("vtk.\n");
      RSS( ref_export_vtk( ref_grid, "validate.vtk" ), "vtk" );

      printf("tec.\n");
      RSS( ref_export_tec( ref_grid, "validate.tec" ), "tec" );

      printf("done.\n");
      return 0;
    }

  return 0;
}
