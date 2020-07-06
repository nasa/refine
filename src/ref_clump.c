
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

#include "ref_clump.h"

#include <stdio.h>
#include <stdlib.h>

#include "ref_adapt.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_malloc.h"
#include "ref_matrix.h"

static REF_STATUS ref_clump_zone_around(FILE *f, REF_CELL ref_cell,
                                        REF_DICT ref_dict, const char *zonetype,
                                        REF_DICT node_dict, REF_NODE ref_node,
                                        REF_INT node) {
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL xyz_comp[3], xyz_phys[3];
  REF_INT local;
  REF_DBL m[6], jacob[9];

  if (ref_dict_n(ref_dict) <= 0) return REF_SUCCESS;

  RSS(ref_node_metric_get(ref_node, node, m), "get");
  RSS(ref_matrix_jacob_m(m, jacob), "jac");

  fprintf(
      f, "zone t=\"%s\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      zonetype, ref_dict_n(node_dict), ref_dict_n(ref_dict), "point", zonetype);

  for (item = 0; item < ref_dict_n(node_dict); item++) {
    local = ref_dict_key(node_dict, item);
    xyz_phys[0] =
        ref_node_xyz(ref_node, 0, local) - ref_node_xyz(ref_node, 0, node);
    xyz_phys[1] =
        ref_node_xyz(ref_node, 1, local) - ref_node_xyz(ref_node, 1, node);
    xyz_phys[2] =
        ref_node_xyz(ref_node, 2, local) - ref_node_xyz(ref_node, 2, node);
    RSS(ref_matrix_vect_mult(jacob, xyz_phys, xyz_comp), "ax");
    fprintf(f, " %.16e %.16e %.16e %.16e %.16e %.16e " REF_GLOB_FMT "\n",
            ref_node_xyz(ref_node, 0, local), ref_node_xyz(ref_node, 1, local),
            ref_node_xyz(ref_node, 2, local), xyz_comp[0], xyz_comp[1],
            xyz_comp[2], ref_node_global(ref_node, local));
  }

  for (item = 0; item < ref_dict_n(ref_dict); item++) {
    cell = ref_dict_key(ref_dict, item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      RSS(ref_dict_location(node_dict, nodes[cell_node], &local), "ret");
      fprintf(f, " %d", local + 1);
    }
    fprintf(f, "\n");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_clump_cell_zone(FILE *f, REF_CELL ref_cell,
                                      REF_DICT ref_dict, const char *zonetype,
                                      REF_NODE ref_node) {
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT local;
  REF_DICT node_dict;

  if (ref_dict_n(ref_dict) <= 0) return REF_SUCCESS;

  RSS(ref_dict_create(&node_dict), "create cell dict");
  for (item = 0; item < ref_dict_n(ref_dict); item++) {
    cell = ref_dict_key(ref_dict, item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
    }
  }

  fprintf(
      f, "zone t=\"%s\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      zonetype, ref_dict_n(node_dict), ref_dict_n(ref_dict), "point", zonetype);

  for (item = 0; item < ref_dict_n(node_dict); item++) {
    local = ref_dict_key(node_dict, item);
    fprintf(f, " %.16e %.16e %.16e\n", ref_node_xyz(ref_node, 0, local),
            ref_node_xyz(ref_node, 1, local), ref_node_xyz(ref_node, 2, local));
  }

  for (item = 0; item < ref_dict_n(ref_dict); item++) {
    cell = ref_dict_key(ref_dict, item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      RSS(ref_dict_location(node_dict, nodes[cell_node], &local), "ret");
      fprintf(f, " %d", local + 1);
    }
    fprintf(f, "\n");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_clump_around(REF_GRID ref_grid, REF_INT node,
                            const char *filename) {
  REF_DICT node_dict, tri_dict, tet_dict;
  REF_DICT ref_dict;
  REF_CELL ref_cell;
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  const char *zonetype;

  FILE *f;

  RSS(ref_dict_create(&node_dict), "create nodes");
  RSS(ref_dict_create(&tri_dict), "create tris");
  RSS(ref_dict_create(&tet_dict), "create tets");

  ref_cell = ref_grid_tri(ref_grid);
  ref_dict = tri_dict;
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    RSS(ref_dict_store(ref_dict, cell, 0), "store");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }
  ref_cell = ref_grid_tet(ref_grid);
  ref_dict = tet_dict;
  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    RSS(ref_dict_store(ref_dict, cell, 0), "store");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }

  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  fprintf(f, "title=\"tecplot refine clump file\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\" \"xm\" \"ym\" \"zm\" \"g\"\n");

  ref_cell = ref_grid_tri(ref_grid);
  ref_dict = tri_dict;
  zonetype = "fetriangle";
  RSS(ref_clump_zone_around(f, ref_cell, ref_dict, zonetype, node_dict,
                            ref_grid_node(ref_grid), node),
      "zone");

  ref_cell = ref_grid_tet(ref_grid);
  ref_dict = tet_dict;
  zonetype = "fetetrahedron";
  RSS(ref_clump_zone_around(f, ref_cell, ref_dict, zonetype, node_dict,
                            ref_grid_node(ref_grid), node),
      "zone");

  fclose(f);

  RSS(ref_dict_free(tet_dict), "free tet");
  RSS(ref_dict_free(tri_dict), "free tris");
  RSS(ref_dict_free(node_dict), "free nodes");

  return REF_SUCCESS;
}

REF_STATUS ref_clump_between(REF_GRID ref_grid, REF_INT node0, REF_INT node1,
                             const char *filename) {
  REF_DICT node_dict, tri_dict, tet_dict;
  REF_DICT ref_dict;
  REF_CELL ref_cell;
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  const char *zonetype;

  FILE *f;

  RSS(ref_dict_create(&node_dict), "create nodes");
  RSS(ref_dict_create(&tri_dict), "create tris");
  RSS(ref_dict_create(&tet_dict), "create tets");

  ref_cell = ref_grid_tri(ref_grid);
  ref_dict = tri_dict;
  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    RSS(ref_dict_store(ref_dict, cell, 0), "store");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }
  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    RSS(ref_dict_store(ref_dict, cell, 0), "store");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }
  ref_cell = ref_grid_tet(ref_grid);
  ref_dict = tet_dict;
  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    RSS(ref_dict_store(ref_dict, cell, 0), "store");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }
  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    RSS(ref_dict_store(ref_dict, cell, 0), "store");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }

  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  fprintf(f, "title=\"tecplot refine clump file\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\" \"xm\" \"ym\" \"zm\" \"g\"\n");

  ref_cell = ref_grid_tri(ref_grid);
  ref_dict = tri_dict;
  zonetype = "fetriangle";
  RSS(ref_clump_zone_around(f, ref_cell, ref_dict, zonetype, node_dict,
                            ref_grid_node(ref_grid), node0),
      "zone");

  ref_cell = ref_grid_tet(ref_grid);
  ref_dict = tet_dict;
  zonetype = "fetetrahedron";
  RSS(ref_clump_zone_around(f, ref_cell, ref_dict, zonetype, node_dict,
                            ref_grid_node(ref_grid), node0),
      "zone");

  fclose(f);

  RSS(ref_dict_free(tet_dict), "free tet");
  RSS(ref_dict_free(tri_dict), "free tris");
  RSS(ref_dict_free(node_dict), "free nodes");

  return REF_SUCCESS;
}

REF_STATUS ref_clump_between_export_to(REF_GRID ref_grid, REF_INT node0,
                                       REF_INT node1, const char *filename) {
  REF_GRID loc_grid;
  REF_NODE ref_node, loc_node;
  REF_CELL ref_cell, loc_cell;
  REF_INT item, cell, cell_node;
  REF_INT old, index, new_cell, new_node, ixyz;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DICT node_dict;

  RSS(ref_grid_create(&loc_grid, ref_grid_mpi(ref_grid)), "create grid");

  RSS(ref_dict_create(&node_dict), "create nodes");

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }
  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }
  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }

  ref_node = ref_grid_node(ref_grid);
  loc_node = ref_grid_node(loc_grid);
  each_ref_dict_key(node_dict, index, old) {
    RSS(ref_node_add(loc_node, index, &new_node), "new_node");
    for (ixyz = 0; ixyz < 3; ixyz++) {
      ref_node_xyz(loc_node, ixyz, new_node) =
          ref_node_xyz(ref_node, ixyz, old);
    }
  }
  RSS(ref_node_initialize_n_global(ref_node, ref_dict_n(node_dict)),
      "init glob");

  ref_cell = ref_grid_tet(ref_grid);
  loc_cell = ref_grid_tet(loc_grid);
  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      old = nodes[cell_node];
      RSS(ref_dict_location(node_dict, old, &(nodes[cell_node])), "map");
    }
    RXS(ref_cell_with(loc_cell, nodes, &new_cell), REF_NOT_FOUND, "exists?");
    if (REF_EMPTY != new_cell) {
      RSS(ref_cell_add(loc_cell, nodes, &new_cell), "new cell");
    }
  }
  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      old = nodes[cell_node];
      RSS(ref_dict_location(node_dict, old, &(nodes[cell_node])), "map");
    }
    RXS(ref_cell_with(loc_cell, nodes, &new_cell), REF_NOT_FOUND, "exists?");
    if (REF_EMPTY != new_cell) {
      RSS(ref_cell_add(loc_cell, nodes, &new_cell), "new cell");
    }
  }

  ref_cell = ref_grid_tri(ref_grid);
  loc_cell = ref_grid_tri(loc_grid);
  each_ref_cell_having_node(ref_cell, node0, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      old = nodes[cell_node];
      RSS(ref_dict_location(node_dict, old, &(nodes[cell_node])), "map");
    }
    RXS(ref_cell_with(loc_cell, nodes, &new_cell), REF_NOT_FOUND, "exists?");
    if (REF_EMPTY != new_cell) {
      RSS(ref_cell_add(loc_cell, nodes, &new_cell), "new cell");
    }
  }
  each_ref_cell_having_node(ref_cell, node1, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      old = nodes[cell_node];
      RSS(ref_dict_location(node_dict, old, &(nodes[cell_node])), "map");
    }
    RXS(ref_cell_with(loc_cell, nodes, &new_cell), REF_NOT_FOUND, "exists?");
    if (REF_EMPTY != new_cell) {
      RSS(ref_cell_add(loc_cell, nodes, &new_cell), "new cell");
    }
  }

  RSS(ref_dict_free(node_dict), "free nodes");

  RSS(ref_export_by_extension(loc_grid, filename), "dump");

  RSS(ref_grid_free(loc_grid), "free loc grid");

  return REF_SUCCESS;
}

REF_STATUS ref_clump_tri_around(REF_GRID ref_grid, REF_INT node,
                                const char *filename) {
  REF_DICT node_dict, tri_dict;
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT local;
  REF_DBL xyz_comp[3], xyz_phys[3];
  REF_DBL m[6], jacob[9];

  FILE *f;

  RSS(ref_dict_create(&node_dict), "create nodes");
  RSS(ref_dict_create(&tri_dict), "create tris");

  each_ref_cell_having_node(ref_cell, node, item, cell) {
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    RSS(ref_dict_store(tri_dict, cell, 0), "store");
    each_ref_cell_cell_node(ref_cell, cell_node)
        RSS(ref_dict_store(node_dict, nodes[cell_node], 0), "store");
  }

  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  fprintf(f, "title=\"tecplot refine clump file\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\"\n");

  fprintf(
      f,
      "zone t=\"clump\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      ref_dict_n(node_dict), ref_dict_n(tri_dict), "point", "fequadrilateral");

  RSS(ref_node_metric_get(ref_grid_node(ref_grid), node, m), "get");
  RSS(ref_matrix_jacob_m(m, jacob), "jac");
  for (item = 0; item < ref_dict_n(node_dict); item++) {
    local = ref_dict_key(node_dict, item);
    xyz_phys[0] = ref_node_xyz(ref_grid_node(ref_grid), 0, local) -
                  ref_node_xyz(ref_grid_node(ref_grid), 0, node);
    xyz_phys[1] = ref_node_xyz(ref_grid_node(ref_grid), 1, local) -
                  ref_node_xyz(ref_grid_node(ref_grid), 1, node);
    xyz_phys[2] = ref_node_xyz(ref_grid_node(ref_grid), 2, local) -
                  ref_node_xyz(ref_grid_node(ref_grid), 2, node);
    RSS(ref_matrix_vect_mult(jacob, xyz_phys, xyz_comp), "ax");
    fprintf(f, " %.16e %.16e %.16e\n", xyz_comp[0], xyz_comp[1], xyz_comp[2]);
  }

  for (item = 0; item < ref_dict_n(tri_dict); item++) {
    cell = ref_dict_key(tri_dict, item);
    RSS(ref_cell_nodes(ref_cell, cell, nodes), "n");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      RSS(ref_dict_location(node_dict, nodes[cell_node], &local), "ret");
      fprintf(f, " %d", local + 1);
    }
    RSS(ref_dict_location(node_dict, nodes[0], &local), "ret");
    fprintf(f, " %d", local + 1);
    fprintf(f, "\n");
  }

  fclose(f);

  RSS(ref_dict_free(tri_dict), "free tris");
  RSS(ref_dict_free(node_dict), "free nodes");

  return REF_SUCCESS;
}

REF_STATUS ref_clump_short_edges(REF_GRID ref_grid, REF_DBL ratio_tol) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_INT ntarget;
  REF_INT node0, node1;
  REF_INT edge;
  REF_DBL edge_ratio;

  char filename[1024];

  RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");

  ntarget = 0;
  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    if (ntarget >= 100) {
      printf("over 100 short edges, giving up\n");
      break;
    }
    node0 = ref_edge_e2n(ref_edge, 0, edge);
    node1 = ref_edge_e2n(ref_edge, 1, edge);
    RSS(ref_node_ratio(ref_node, node0, node1, &edge_ratio), "ratio");
    if (edge_ratio < ratio_tol * ref_grid_adapt(ref_grid, collapse_ratio)) {
      sprintf(filename, "clump%d.t", ntarget);
      RSS(ref_clump_between(ref_grid, node0, node1, filename), "dump");
      ntarget++;
    }
  }

  RSS(ref_edge_free(ref_edge), "free edges");

  return REF_SUCCESS;
}

REF_STATUS ref_clump_short_edges_twod(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT ntarget;
  REF_INT node, node0, node1;
  REF_INT edge;
  REF_DBL edge_ratio;

  char filename[1024];

  RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");

  ref_malloc_init(ratio, ref_node_max(ref_node), REF_DBL,
                  2.0 * ref_grid_adapt(ref_grid, collapse_ratio));

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    node0 = ref_edge_e2n(ref_edge, 0, edge);
    node1 = ref_edge_e2n(ref_edge, 1, edge);

    RSS(ref_node_ratio(ref_node, node0, node1, &edge_ratio), "ratio");
    ratio[node0] = MIN(ratio[node0], edge_ratio);
    ratio[node1] = MIN(ratio[node1], edge_ratio);
  }

  ntarget = 0;
  for (node = 0; node < ref_node_max(ref_node); node++)
    if (ratio[node] < ref_grid_adapt(ref_grid, collapse_ratio)) {
      sprintf(filename, "clump%d.t", ntarget);
      RSS(ref_clump_tri_around(ref_grid, node, filename), "dump");
      ntarget++;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_clump_long_edges(REF_GRID ref_grid, REF_DBL ratio_tol) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_INT ntarget;
  REF_INT node0, node1;
  REF_INT edge;
  REF_DBL edge_ratio;

  char filename[1024];

  RSS(ref_edge_create(&ref_edge, ref_grid), "orig edges");

  ntarget = 0;
  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    if (ntarget >= 100) {
      printf("over 100 short edges, giving up\n");
      break;
    }
    node0 = ref_edge_e2n(ref_edge, 0, edge);
    node1 = ref_edge_e2n(ref_edge, 1, edge);
    RSS(ref_node_ratio(ref_node, node0, node1, &edge_ratio), "ratio");
    if (edge_ratio > ratio_tol * ref_grid_adapt(ref_grid, post_max_ratio)) {
      sprintf(filename, "clump%d.t", ntarget);
      RSS(ref_clump_between(ref_grid, node0, node1, filename), "dump");
      ntarget++;
    }
  }

  RSS(ref_edge_free(ref_edge), "free edges");

  return REF_SUCCESS;
}

REF_STATUS ref_clump_tet_quality(REF_GRID ref_grid, REF_DBL min_quality,
                                 const char *filename) {
  REF_DICT ref_dict;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL quality;
  const char *zonetype;
  FILE *f;

  RSS(ref_dict_create(&ref_dict), "create cell dict");

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    if (quality < min_quality) RSS(ref_dict_store(ref_dict, cell, 0), "store");
  }

  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  fprintf(f, "title=\"tecplot refine cell clump file\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\"\n");

  zonetype = "fetetrahedron";
  RSS(ref_clump_cell_zone(f, ref_cell, ref_dict, zonetype, ref_node), "zone");

  fclose(f);

  RSS(ref_dict_free(ref_dict), "free tet");

  return REF_SUCCESS;
}
