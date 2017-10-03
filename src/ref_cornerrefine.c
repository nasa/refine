
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

#include "ref_inflate.h"

#include "ref_math.h"

#include "ref_import.h"
#include "ref_export.h"

#include "ref_cell.h"
#include "ref_grid.h"
#include "ref_sort.h"
#include "ref_adj.h"
#include "ref_node.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_dict.h"
#include "ref_list.h"

#include "ref_edge.h"

#include "ref_fixture.h"
#include "ref_metric.h"
#include "ref_gather.h"

#include "ref_subdiv.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_CELL ref_cell;
  REF_NODE ref_node;
  REF_SUBDIV ref_subdiv;
  REF_BOOL valid_inputs;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node, item;
  REF_INT arg_parse_position;
  REF_INT best_node;
  REF_DBL x_node, y_node, best_dist;
  REF_DBL dx, dy, dist;

  valid_inputs = ( 3 <= argc );

  if ( !valid_inputs )
    {
      printf("usage: %s input_grid.extension output_grid.extension -xy x_node y_node ...\n",argv[0]);
      return 1;
    }

  RSS( ref_import_by_extension( &ref_grid, argv[1] ), "in");
  RSS( ref_metric_unit_node( ref_grid_node(ref_grid) ), "met");

  ref_cell = ref_grid_pri(ref_grid);
  ref_node = ref_grid_node(ref_grid);

  RSS( ref_subdiv_create( &ref_subdiv, ref_grid ), "init" );
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
  {
    RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[1], nodes[2] ), "o0" );
    RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[2], nodes[0] ), "o1" );
    RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[0], nodes[1] ), "o2" );
  }
  RSS(ref_subdiv_split(ref_subdiv),"split");
  RSS(ref_subdiv_free(ref_subdiv),"free");

  arg_parse_position = 3;
  while ( arg_parse_position < argc &&
          strcmp(argv[arg_parse_position],"-xy") == 0 )
    {
      if (argc < arg_parse_position+3)
        THROW("incomplete -xy argument");
      printf("%s %s %s\n",
             argv[arg_parse_position],
             argv[arg_parse_position+1],
             argv[arg_parse_position+2]);
      x_node = atof(argv[arg_parse_position+1]);
      y_node = atof(argv[arg_parse_position+2]);
      printf("%s %.15e %.15e\n",
             argv[arg_parse_position],
             x_node,
             y_node);
      arg_parse_position += 3;
      best_dist = 0;
      best_node = REF_EMPTY;
      each_ref_node_valid_node(ref_node,node)
      {
        dx = ref_node_xyz(ref_node,0,node)-x_node;
        dy = ref_node_xyz(ref_node,2,node)-y_node;
        dist = sqrt(dx*dx+dy*dy);
        if ( best_node == REF_EMPTY )
          {
            best_node = node;
            best_dist = dist;
          }
        else
          {
            if ( dist < best_dist )
              {
                best_node = node;
                best_dist = dist;
              }
          }
      }

      if ( best_node == REF_EMPTY )
        THROW("missing nodes?");

      printf("node %d x %.15e y %.15e dist %.15e\n",
             best_node,
             ref_node_xyz(ref_node,0,best_node),
             ref_node_xyz(ref_node,2,best_node),
             best_dist);

      RSS( ref_subdiv_create( &ref_subdiv, ref_grid ), "init" );
      node = best_node;
      each_ref_cell_having_node( ref_cell, node, item, cell )
      {
        RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
        RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[1], nodes[2] ), "o0" );
        RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[2], nodes[0] ), "o1" );
        RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[0], nodes[1] ), "o2" );
      }
      RSS(ref_subdiv_split(ref_subdiv),"split");
      RSS(ref_subdiv_free(ref_subdiv),"free");
    }

  RSS( ref_export_by_extension( ref_grid, argv[2] ), "out");

  RSS(ref_grid_free(ref_grid),"free");

  return 0;
}
