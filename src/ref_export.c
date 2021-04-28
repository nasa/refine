
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

#include "ref_export.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_endian.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"

#define VTK_TRIANGLE (5)
#define VTK_QUAD (9)
#define VTK_TETRA (10)
#define VTK_HEXAHEDRON (12)
#define VTK_WEDGE (13)
#define VTK_PYRAMID (14)

/*
3-4 UGRID
| |\
| | 2
| |/
0-1
2-3 VTK
| |\
| | 4
| |/
1-0
 */

#define VTK_PYRAMID_ORDER(vtk_nodes) \
  {                                  \
    REF_INT ugrid_nodes[5];          \
    ugrid_nodes[0] = (vtk_nodes)[0]; \
    ugrid_nodes[1] = (vtk_nodes)[1]; \
    ugrid_nodes[2] = (vtk_nodes)[2]; \
    ugrid_nodes[3] = (vtk_nodes)[3]; \
    ugrid_nodes[4] = (vtk_nodes)[4]; \
    (vtk_nodes)[0] = ugrid_nodes[1]; \
    (vtk_nodes)[1] = ugrid_nodes[0]; \
    (vtk_nodes)[2] = ugrid_nodes[3]; \
    (vtk_nodes)[3] = ugrid_nodes[4]; \
    (vtk_nodes)[4] = ugrid_nodes[2]; \
  }

/*
 /3-/0
5-+-2| UGRID
 \4-\1
 /4-/1
5-+-2| VTK
 \3-\0
*/
#define VTK_WEDGE_ORDER(vtk_nodes)   \
  {                                  \
    REF_INT ugrid_nodes[6];          \
    ugrid_nodes[0] = (vtk_nodes)[0]; \
    ugrid_nodes[1] = (vtk_nodes)[1]; \
    ugrid_nodes[2] = (vtk_nodes)[2]; \
    ugrid_nodes[3] = (vtk_nodes)[3]; \
    ugrid_nodes[4] = (vtk_nodes)[4]; \
    ugrid_nodes[5] = (vtk_nodes)[5]; \
    (vtk_nodes)[0] = ugrid_nodes[1]; \
    (vtk_nodes)[1] = ugrid_nodes[0]; \
    (vtk_nodes)[2] = ugrid_nodes[2]; \
    (vtk_nodes)[3] = ugrid_nodes[4]; \
    (vtk_nodes)[4] = ugrid_nodes[3]; \
    (vtk_nodes)[5] = ugrid_nodes[5]; \
  }

static REF_STATUS ref_export_faceid_range(REF_GRID ref_grid,
                                          REF_INT *min_faceid,
                                          REF_INT *max_faceid) {
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  *min_faceid = REF_INT_MAX;
  *max_faceid = REF_INT_MIN;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    *min_faceid = MIN(*min_faceid, nodes[ref_cell_id_index(ref_cell)]);
    *max_faceid = MAX(*max_faceid, nodes[ref_cell_id_index(ref_cell)]);
  }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    *min_faceid = MIN(*min_faceid, nodes[ref_cell_id_index(ref_cell)]);
    *max_faceid = MAX(*max_faceid, nodes[ref_cell_id_index(ref_cell)]);
  }

  return REF_SUCCESS;
}
static REF_STATUS ref_export_cell_id_range(REF_CELL ref_cell,
                                           REF_INT *min_faceid,
                                           REF_INT *max_faceid) {
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  *min_faceid = REF_INT_MAX;
  *max_faceid = REF_INT_MIN;

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    *min_faceid = MIN(*min_faceid, nodes[ref_cell_id_index(ref_cell)]);
    *max_faceid = MAX(*max_faceid, nodes[ref_cell_id_index(ref_cell)]);
  }

  return REF_SUCCESS;
}

/* https://www.vtk.org/VTK/img/file-formats.pdf */
static REF_STATUS ref_export_vtk(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT ncell, size;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell;
  REF_INT group;

  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "# vtk DataFile Version 2.0\n");
  fprintf(file, "ref_export_vtk\n");
  fprintf(file, "ASCII\n");

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  fprintf(file, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(file, "POINTS %d double\n", ref_node_n(ref_node));

  for (node = 0; node < ref_node_n(ref_node); node++) {
    fprintf(file, " %.16e %.16e %.16e\n", ref_node_xyz(ref_node, 0, n2o[node]),
            ref_node_xyz(ref_node, 1, n2o[node]),
            ref_node_xyz(ref_node, 2, n2o[node]));
  }

  ncell = 0;
  size = 0;

  each_ref_grid_2d_3d_ref_cell(ref_grid, group, ref_cell) {
    ncell += ref_cell_n(ref_cell);
    size += ref_cell_n(ref_cell) * (1 + ref_cell_node_per(ref_cell));
  }

  fprintf(file, "CELLS %d %d\n", ncell, size);

  each_ref_grid_2d_3d_ref_cell(ref_grid, group, ref_cell) {
    node_per = ref_cell_node_per(ref_cell);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      fprintf(file, " %d", node_per);
      if (5 == node_per) VTK_PYRAMID_ORDER(nodes);
      if (6 == node_per) VTK_WEDGE_ORDER(nodes);
      for (node = 0; node < node_per; node++)
        fprintf(file, " %d", o2n[nodes[node]]);
      fprintf(file, "\n");
    }
  }

  fprintf(file, "CELL_TYPES %d\n", ncell);

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) fprintf(file, " %d\n", VTK_TRIANGLE);

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) fprintf(file, " %d\n", VTK_QUAD);

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) fprintf(file, " %d\n", VTK_TETRA);

  ref_cell = ref_grid_pyr(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) fprintf(file, " %d\n", VTK_PYRAMID);

  ref_cell = ref_grid_pri(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) fprintf(file, " %d\n", VTK_WEDGE);

  ref_cell = ref_grid_hex(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell)
      fprintf(file, " %d\n", VTK_HEXAHEDRON);

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_export_tec_edg_zone(REF_GRID ref_grid, FILE *file) {
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *g2l, *l2g;
  REF_INT nedge;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell, cell_node;
  REF_INT nnode;
  REF_DICT ref_dict;
  REF_INT boundary_tag, boundary_index;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_dict_create(&ref_dict), "create dict");

  ref_cell = ref_grid_edg(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_dict_store(ref_dict, nodes[ref_cell_id_index(ref_cell)], REF_EMPTY),
        "mark tri");
  }

  each_ref_dict_key(ref_dict, boundary_index, boundary_tag) {
    RSS(ref_grid_cell_id_nodes(ref_grid, ref_cell, boundary_tag, &nnode, &nedge,
                               &g2l, &l2g),
        "extract this edge");

    fprintf(file,
            "zone t=\"edge%d\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
            boundary_tag, nnode, nedge, "point", "felineseg");

    for (node = 0; node < nnode; node++)
      fprintf(file, " %.16e %.16e %.16e\n",
              ref_node_xyz(ref_node, 0, l2g[node]),
              ref_node_xyz(ref_node, 1, l2g[node]),
              ref_node_xyz(ref_node, 2, l2g[node]));

    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (boundary_tag == nodes[ref_cell_id_index(ref_cell)]) {
        each_ref_cell_cell_node(ref_cell, cell_node) {
          fprintf(file, " %d", g2l[nodes[cell_node]] + 1);
        }
        fprintf(file, "\n");
      }
    }

    ref_free(l2g);
    ref_free(g2l);
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

static REF_STATUS ref_export_tec_ed2_zone(REF_GRID ref_grid, FILE *file) {
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *g2l, *l2g;
  REF_INT nedge;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnode;
  REF_DICT ref_dict;
  REF_INT boundary_tag, boundary_index;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_dict_create(&ref_dict), "create dict");

  ref_cell = ref_grid_ed2(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_dict_store(ref_dict, nodes[ref_cell_id_index(ref_cell)], REF_EMPTY),
        "mark tri");
  }

  each_ref_dict_key(ref_dict, boundary_index, boundary_tag) {
    RSS(ref_grid_cell_id_nodes(ref_grid, ref_cell, boundary_tag, &nnode, &nedge,
                               &g2l, &l2g),
        "extract this edge");

    fprintf(file,
            "zone t=\"p2edge%d\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
            boundary_tag, nnode, 2 * nedge, "point", "felineseg");

    for (node = 0; node < nnode; node++)
      fprintf(file, " %.16e %.16e %.16e\n",
              ref_node_xyz(ref_node, 0, l2g[node]),
              ref_node_xyz(ref_node, 1, l2g[node]),
              ref_node_xyz(ref_node, 2, l2g[node]));

    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (boundary_tag == nodes[ref_cell_id_index(ref_cell)]) {
        fprintf(file, " %d %d\n", g2l[nodes[0]] + 1, g2l[nodes[2]] + 1);
        fprintf(file, " %d %d\n", g2l[nodes[2]] + 1, g2l[nodes[1]] + 1);
      }
    }

    ref_free(l2g);
    ref_free(g2l);
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

static REF_STATUS ref_export_tec_ed3_zone(REF_GRID ref_grid, FILE *file) {
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *g2l, *l2g;
  REF_INT nedge;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnode;
  REF_DICT ref_dict;
  REF_INT boundary_tag, boundary_index;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_dict_create(&ref_dict), "create dict");

  ref_cell = ref_grid_ed3(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_dict_store(ref_dict, nodes[ref_cell_id_index(ref_cell)], REF_EMPTY),
        "mark tri");
  }

  each_ref_dict_key(ref_dict, boundary_index, boundary_tag) {
    RSS(ref_grid_cell_id_nodes(ref_grid, ref_cell, boundary_tag, &nnode, &nedge,
                               &g2l, &l2g),
        "extract this edge");

    fprintf(file,
            "zone t=\"p3edge%d\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
            boundary_tag, nnode, 3 * nedge, "point", "felineseg");

    for (node = 0; node < nnode; node++)
      fprintf(file, " %.16e %.16e %.16e\n",
              ref_node_xyz(ref_node, 0, l2g[node]),
              ref_node_xyz(ref_node, 1, l2g[node]),
              ref_node_xyz(ref_node, 2, l2g[node]));

    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (boundary_tag == nodes[ref_cell_id_index(ref_cell)]) {
        fprintf(file, " %d %d\n", g2l[nodes[0]] + 1, g2l[nodes[2]] + 1);
        fprintf(file, " %d %d\n", g2l[nodes[2]] + 1, g2l[nodes[3]] + 1);
        fprintf(file, " %d %d\n", g2l[nodes[3]] + 1, g2l[nodes[1]] + 1);
      }
    }

    ref_free(l2g);
    ref_free(g2l);
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

static REF_STATUS ref_export_tec_surf_zone(REF_GRID ref_grid, FILE *file) {
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *g2l, *l2g;
  REF_INT nface;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnode;
  REF_DICT ref_dict;
  REF_INT boundary_tag, boundary_index;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_dict_create(&ref_dict), "create dict");

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) RSS(
      ref_dict_store(ref_dict, nodes[ref_cell_id_index(ref_cell)], REF_EMPTY),
      "mark tri");

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) RSS(
      ref_dict_store(ref_dict, nodes[ref_cell_id_index(ref_cell)], REF_EMPTY),
      "mark qua");

  each_ref_dict_key(ref_dict, boundary_index, boundary_tag) {
    RSS(ref_grid_tri_qua_id_nodes(ref_grid, boundary_tag, &nnode, &nface, &g2l,
                                  &l2g),
        "extract this boundary");

    fprintf(file,
            "zone t=\"face%d\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
            boundary_tag, nnode, nface, "point", "fequadrilateral");

    for (node = 0; node < nnode; node++)
      fprintf(file, " %.16e %.16e %.16e\n",
              ref_node_xyz(ref_node, 0, l2g[node]),
              ref_node_xyz(ref_node, 1, l2g[node]),
              ref_node_xyz(ref_node, 2, l2g[node]));

    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (boundary_tag == nodes[ref_cell_id_index(ref_cell)]) {
        nodes[3] = nodes[2];
        for (node = 0; node < 4; node++) {
          fprintf(file, " %d", g2l[nodes[node]] + 1);
        }
        fprintf(file, "\n");
      }
    }

    ref_cell = ref_grid_qua(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (boundary_tag == nodes[ref_cell_id_index(ref_cell)]) {
        for (node = 0; node < 4; node++)
          fprintf(file, " %d", g2l[nodes[node]] + 1);
        fprintf(file, "\n");
      }
    }

    ref_free(l2g);
    ref_free(g2l);
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

static REF_STATUS ref_export_tec_tr2_zone(REF_GRID ref_grid, FILE *file) {
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *g2l, *l2g;
  REF_INT ncell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnode;
  REF_DICT ref_dict;
  REF_INT boundary_tag, boundary_index;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_dict_create(&ref_dict), "create dict");

  ref_cell = ref_grid_tr2(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) RSS(
      ref_dict_store(ref_dict, nodes[ref_cell_id_index(ref_cell)], REF_EMPTY),
      "mark tri");

  each_ref_dict_key(ref_dict, boundary_index, boundary_tag) {
    RSS(ref_grid_cell_id_nodes(ref_grid, ref_cell, boundary_tag, &nnode, &ncell,
                               &g2l, &l2g),
        "extract this tri");

    fprintf(file,
            "zone t=\"p2face%d\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
            boundary_tag, nnode, 4 * ncell, "point", "fetriangle");

    for (node = 0; node < nnode; node++)
      fprintf(file, " %.16e %.16e %.16e\n",
              ref_node_xyz(ref_node, 0, l2g[node]),
              ref_node_xyz(ref_node, 1, l2g[node]),
              ref_node_xyz(ref_node, 2, l2g[node]));

    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (boundary_tag == nodes[ref_cell_id_index(ref_cell)]) {
        fprintf(file, " %d %d %d\n", g2l[nodes[0]] + 1, g2l[nodes[3]] + 1,
                g2l[nodes[5]] + 1);
        fprintf(file, " %d %d %d\n", g2l[nodes[1]] + 1, g2l[nodes[4]] + 1,
                g2l[nodes[3]] + 1);
        fprintf(file, " %d %d %d\n", g2l[nodes[2]] + 1, g2l[nodes[5]] + 1,
                g2l[nodes[4]] + 1);
        fprintf(file, " %d %d %d\n", g2l[nodes[3]] + 1, g2l[nodes[4]] + 1,
                g2l[nodes[5]] + 1);
      }
    }

    ref_free(l2g);
    ref_free(g2l);
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

static REF_STATUS ref_export_tec_tr3_zone(REF_GRID ref_grid, FILE *file) {
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *g2l, *l2g;
  REF_INT ncell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnode;
  REF_DICT ref_dict;
  REF_INT boundary_tag, boundary_index;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_dict_create(&ref_dict), "create dict");

  ref_cell = ref_grid_tr3(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) RSS(
      ref_dict_store(ref_dict, nodes[ref_cell_id_index(ref_cell)], REF_EMPTY),
      "mark tri");

  each_ref_dict_key(ref_dict, boundary_index, boundary_tag) {
    RSS(ref_grid_cell_id_nodes(ref_grid, ref_cell, boundary_tag, &nnode, &ncell,
                               &g2l, &l2g),
        "extract this tri");

    fprintf(file,
            "zone t=\"p3face%d\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
            boundary_tag, nnode, 9 * ncell, "point", "fetriangle");

    for (node = 0; node < nnode; node++)
      fprintf(file, " %.16e %.16e %.16e\n",
              ref_node_xyz(ref_node, 0, l2g[node]),
              ref_node_xyz(ref_node, 1, l2g[node]),
              ref_node_xyz(ref_node, 2, l2g[node]));

    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (boundary_tag == nodes[ref_cell_id_index(ref_cell)]) {
        fprintf(file, " %d %d %d\n", g2l[nodes[0]] + 1, g2l[nodes[3]] + 1,
                g2l[nodes[8]] + 1);
        fprintf(file, " %d %d %d\n", g2l[nodes[8]] + 1, g2l[nodes[3]] + 1,
                g2l[nodes[9]] + 1);

        fprintf(file, " %d %d %d\n", g2l[nodes[1]] + 1, g2l[nodes[5]] + 1,
                g2l[nodes[4]] + 1);
        fprintf(file, " %d %d %d\n", g2l[nodes[4]] + 1, g2l[nodes[5]] + 1,
                g2l[nodes[9]] + 1);

        fprintf(file, " %d %d %d\n", g2l[nodes[2]] + 1, g2l[nodes[7]] + 1,
                g2l[nodes[6]] + 1);
        fprintf(file, " %d %d %d\n", g2l[nodes[6]] + 1, g2l[nodes[7]] + 1,
                g2l[nodes[9]] + 1);

        fprintf(file, " %d %d %d\n", g2l[nodes[9]] + 1, g2l[nodes[5]] + 1,
                g2l[nodes[6]] + 1);
        fprintf(file, " %d %d %d\n", g2l[nodes[9]] + 1, g2l[nodes[7]] + 1,
                g2l[nodes[8]] + 1);
        fprintf(file, " %d %d %d\n", g2l[nodes[9]] + 1, g2l[nodes[3]] + 1,
                g2l[nodes[4]] + 1);
      }
    }

    ref_free(l2g);
    ref_free(g2l);
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

static REF_STATUS ref_export_tec_twod_zone(REF_GRID ref_grid, FILE *file) {
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *g2l, *l2g;
  REF_INT nface;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnode;
  REF_DICT ref_dict;
  REF_INT boundary_tag, boundary_index;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_dict_create(&ref_dict), "create dict");

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) RSS(
      ref_dict_store(ref_dict, nodes[ref_cell_id_index(ref_cell)], REF_EMPTY),
      "mark tri");

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) RSS(
      ref_dict_store(ref_dict, nodes[ref_cell_id_index(ref_cell)], REF_EMPTY),
      "mark qua");

  each_ref_dict_key(ref_dict, boundary_index, boundary_tag) {
    RSS(ref_grid_tri_qua_id_nodes(ref_grid, boundary_tag, &nnode, &nface, &g2l,
                                  &l2g),
        "extract this boundary");

    fprintf(file,
            "zone t=\"face%d\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
            boundary_tag, nnode, nface, "point", "fequadrilateral");

    for (node = 0; node < nnode; node++)
      fprintf(file, " %.16e %.16e\n", ref_node_xyz(ref_node, 0, l2g[node]),
              ref_node_xyz(ref_node, 1, l2g[node]));

    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (boundary_tag == nodes[ref_cell_id_index(ref_cell)]) {
        nodes[3] = nodes[2];
        for (node = 0; node < 4; node++) {
          fprintf(file, " %d", g2l[nodes[node]] + 1);
        }
        fprintf(file, "\n");
      }
    }

    ref_cell = ref_grid_qua(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (boundary_tag == nodes[ref_cell_id_index(ref_cell)]) {
        for (node = 0; node < 4; node++)
          fprintf(file, " %d", g2l[nodes[node]] + 1);
        fprintf(file, "\n");
      }
    }

    ref_free(l2g);
    ref_free(g2l);
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

static REF_STATUS ref_export_tec_vol_zone(REF_GRID ref_grid, FILE *file) {
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT brick[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT nnode;
  REF_INT group, node_per;

  ref_node = ref_grid_node(ref_grid);

  each_ref_grid_3d_ref_cell(ref_grid, group, ref_cell) {
    if (ref_cell_n(ref_cell) > 0) {
      node_per = ref_cell_node_per(ref_cell);

      ref_malloc(o2n, ref_node_max(ref_node), REF_INT);

      for (node = 0; node < ref_node_max(ref_node); node++)
        o2n[node] = REF_EMPTY;

      /* mark nodes needed by this element type */
      nnode = 0;
      each_ref_cell_valid_cell_with_nodes(
          ref_cell, cell,
          nodes) for (node = 0; node < node_per;
                      node++) if (REF_EMPTY == o2n[nodes[node]]) {
        o2n[nodes[node]] = nnode;
        nnode++;
      }

      /* retain node numbering if grid is compact and single element type */
      nnode = 0;
      for (node = 0; node < ref_node_max(ref_node); node++)
        if (REF_EMPTY != o2n[node]) {
          o2n[node] = nnode;
          nnode++;
        }

      ref_malloc(n2o, nnode, REF_INT);

      for (node = 0; node < ref_node_max(ref_node); node++)
        if (REF_EMPTY != o2n[node]) n2o[o2n[node]] = node;

      fprintf(file,
              "zone t=\"e%d\", nodes=%d, elements=%d, datapacking=%s, "
              "zonetype=%s\n",
              node_per, nnode, ref_cell_n(ref_cell), "point", "febrick");

      for (node = 0; node < nnode; node++)
        fprintf(file, " %.16e %.16e %.16e\n",
                ref_node_xyz(ref_node, 0, n2o[node]),
                ref_node_xyz(ref_node, 1, n2o[node]),
                ref_node_xyz(ref_node, 2, n2o[node]));

      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        switch (ref_cell_node_per(ref_cell)) {
          case 4:
            REF_CELL_TEC_BRICK_TET(brick, nodes);
            break;
          case 5:
            REF_CELL_TEC_BRICK_PYR(brick, nodes);
            break;
          case 6:
            REF_CELL_TEC_BRICK_PRI(brick, nodes);
            break;
          case 8:
            REF_CELL_TEC_BRICK_HEX(brick, nodes);
            break;
          default:
            RSS(REF_IMPLEMENT, "wrong nodes per cell");
            break;
        }

        for (node = 0; node < 8; node++) {
          fprintf(file, " %d", o2n[brick[node]] + 1);
        }
        fprintf(file, "\n");
      }

      ref_free(n2o);
      ref_free(o2n);
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_export_tec(REF_GRID ref_grid, const char *filename) {
  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine geometry file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

  RSS(ref_export_tec_vol_zone(ref_grid, file), "vol");
  RSS(ref_export_tec_surf_zone(ref_grid, file), "surf");
  RSS(ref_export_tec_tr2_zone(ref_grid, file), "tr2");
  RSS(ref_export_tec_tr3_zone(ref_grid, file), "tr3");
  RSS(ref_export_tec_edg_zone(ref_grid, file), "edg");
  RSS(ref_export_tec_ed2_zone(ref_grid, file), "ed2");
  RSS(ref_export_tec_ed3_zone(ref_grid, file), "ed3");

  fclose(file);
  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_surf(REF_GRID ref_grid, const char *filename) {
  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine geometry file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

  RSS(ref_export_tec_surf_zone(ref_grid, file), "surf");

  fclose(file);
  return REF_SUCCESS;
}

static REF_STATUS ref_export_tec_metric_ellipse_twod(
    REF_GRID ref_grid, const char *root_filename) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT nnode, node;
  REF_INT *o2n, *n2o;
  REF_INT ncell, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL m[6], d[12];
  REF_DBL x, y;
  REF_DBL ex, ey;
  FILE *file;
  char viz_file[256];
  REF_INT i, n = 12;
  REF_INT e0, e1, eb;
  REF_DBL best_z;
  REF_DBL dt = ref_math_in_radians(360.0 / (REF_DBL)n);
  REF_DBL scale = 0.5; /* so the ellipses touch for an ideal grid */

  sprintf(viz_file, "%s_n%d_p%d_ellipse.tec", root_filename,
          ref_mpi_n(ref_grid_mpi(ref_grid)),
          ref_mpi_rank(ref_grid_mpi(ref_grid)));

  file = fopen(viz_file, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", viz_file);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine metric axes\"\n");
  fprintf(file, "variables = \"x\" \"y\"\n");

  ref_malloc_init(o2n, ref_node_max(ref_node), REF_INT, REF_EMPTY);
  nnode = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (REF_EMPTY == o2n[nodes[node]]) {
        o2n[nodes[node]] = nnode;
        nnode++;
      }
    }
  }

  ref_malloc(n2o, nnode, REF_INT);
  for (node = 0; node < ref_node_max(ref_node); node++) {
    if (REF_EMPTY != o2n[node]) {
      n2o[o2n[node]] = node;
    }
  }

  ncell = nnode * n;

  fprintf(file,
          "zone t=\"ellipse\", nodes=%d, elements=%d, datapacking=%s, "
          "zonetype=%s\n",
          1 * ncell, 1 * ncell, "point", "felineseg");

  for (node = 0; node < nnode; node++) {
    RSS(ref_node_metric_get(ref_node, n2o[node], m), "get");
    RSS(ref_matrix_diag_m(m, d), "diag");
    RSS(ref_matrix_ascending_eig(d), "sort eig");
    eb = REF_EMPTY;
    best_z = -1.0;
    for (e0 = 0; e0 < 3; e0++) {
      if (ABS(ref_matrix_vec(d, 2, e0)) > best_z) {
        eb = e0;
        best_z = ABS(ref_matrix_vec(d, 2, e0));
      }
    }
    RUS(REF_EMPTY, eb, "couldn't find best y");
    e0 = eb + 1;
    if (e0 > 2) e0 -= 3;
    e1 = eb + 2;
    if (e1 > 2) e1 -= 3;
    for (i = 0; i < n; i++) {
      ex = scale * cos(i * dt) / sqrt(d[e0]);
      ey = scale * sin(i * dt) / sqrt(d[e1]);
      x = d[3 + 3 * e0] * ex + d[3 + 3 * e1] * ey;
      y = d[4 + 3 * e0] * ex + d[4 + 3 * e1] * ey;
      fprintf(file, " %.16e %.16e\n", ref_node_xyz(ref_node, 0, n2o[node]) + x,
              ref_node_xyz(ref_node, 1, n2o[node]) + y);
    }
  }

  e0 = 0;
  for (node = 0; node < nnode; node++) {
    for (i = 0; i < n - 1; i++)
      fprintf(file, " %d %d\n", i + node * n + 1, i + 1 + node * n + 1);
    fprintf(file, " %d %d\n", n + node * n, 1 + node * n);
  }

  ref_free(n2o);
  ref_free(o2n);

  RSS(ref_export_tec_twod_zone(ref_grid, file), "ellipse surf");

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_metric_ellipse(REF_GRID ref_grid,
                                         const char *root_filename) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT nnode, node;
  REF_INT *o2n, *n2o;
  REF_INT ncell, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL m[6], d[12];
  REF_DBL x, y, z;
  REF_DBL ex, ey;
  FILE *file;
  char viz_file[256];
  REF_INT i, n = 36;
  REF_INT e0, e1;
  REF_DBL dt = ref_math_in_radians(360.0 / (REF_DBL)n);
  REF_DBL scale = 0.5; /* so the ellipses touch for an ideal grid */

  if (ref_grid_twod(ref_grid)) {
    RSS(ref_export_tec_metric_ellipse_twod(ref_grid, root_filename), "2d met");
    return REF_SUCCESS;
  }

  sprintf(viz_file, "%s_n%d_p%d_ellipse.tec", root_filename,
          ref_mpi_n(ref_grid_mpi(ref_grid)),
          ref_mpi_rank(ref_grid_mpi(ref_grid)));

  file = fopen(viz_file, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", viz_file);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine metric axes\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

  ref_malloc_init(o2n, ref_node_max(ref_node), REF_INT, REF_EMPTY);
  nnode = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      if (REF_EMPTY == o2n[nodes[node]]) {
        o2n[nodes[node]] = nnode;
        nnode++;
      }
    }
  }

  ref_malloc(n2o, nnode, REF_INT);
  for (node = 0; node < ref_node_max(ref_node); node++) {
    if (REF_EMPTY != o2n[node]) {
      n2o[o2n[node]] = node;
    }
  }

  ncell = nnode * n;

  fprintf(file,
          "zone t=\"ellipse\", nodes=%d, elements=%d, datapacking=%s, "
          "zonetype=%s\n",
          3 * ncell, 3 * ncell, "point", "felineseg");

  for (node = 0; node < nnode; node++) {
    RSS(ref_node_metric_get(ref_node, n2o[node], m), "get");
    RSS(ref_matrix_diag_m(m, d), "diag");
    RSS(ref_matrix_ascending_eig(d), "sort eig");
    for (e0 = 0; e0 < 3; e0++) {
      e1 = e0 + 1;
      if (e1 == 3) e1 = 0;
      for (i = 0; i < n; i++) {
        ex = scale * cos(i * dt) / sqrt(d[e0]);
        ey = scale * sin(i * dt) / sqrt(d[e1]);
        x = d[3 + 3 * e0] * ex + d[3 + 3 * e1] * ey;
        y = d[4 + 3 * e0] * ex + d[4 + 3 * e1] * ey;
        z = d[5 + 3 * e0] * ex + d[5 + 3 * e1] * ey;
        fprintf(file, " %.16e %.16e %.16e\n",
                ref_node_xyz(ref_node, 0, n2o[node]) + x,
                ref_node_xyz(ref_node, 1, n2o[node]) + y,
                ref_node_xyz(ref_node, 2, n2o[node]) + z);
      }
    }
  }

  for (e0 = 0; e0 < 3; e0++)
    for (node = 0; node < nnode; node++) {
      for (i = 0; i < n - 1; i++)
        fprintf(file, " %d %d\n", i + node * n + 1 + ncell * e0,
                i + 1 + node * n + 1 + ncell * e0);
      fprintf(file, " %d %d\n", n + node * n + ncell * e0,
              1 + node * n + ncell * e0);
    }

  ref_free(n2o);
  ref_free(o2n);

  RSS(ref_export_tec_surf_zone(ref_grid, file), "ellipse surf");

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_export_poly(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nnode, ntri;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell;
  REF_INT dim, attr, mark;
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  /* Part 1 - node list */
  nnode = ref_node_n(ref_node);
  dim = 3;
  attr = 0;
  mark = 0;
  fprintf(file, "%d  %d  %d  %d\n", nnode, dim, attr, mark);

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  for (node = 0; node < ref_node_n(ref_node); node++)
    fprintf(file, "%d  %.16e  %.16e  %.16e\n", node,
            ref_node_xyz(ref_node, 0, n2o[node]),
            ref_node_xyz(ref_node, 1, n2o[node]),
            ref_node_xyz(ref_node, 2, n2o[node]));

  /* Part 2 - facet list */
  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  ntri = ref_cell_n(ref_cell);
  mark = 1;
  fprintf(file, "%d  %d\n", ntri, mark);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    fprintf(file, "1 0 1\n%d  %d %d %d  %d\n", node_per, o2n[nodes[0]],
            o2n[nodes[1]], o2n[nodes[2]], nodes[3]);
  }

  /* Part 3 hole - list */
  /* create one hole point for each convex edge,
     expect right hand triangle to point into domain */
  {
    REF_INT faceid, min_faceid, max_faceid;
    REF_INT nhole, largest_triangle;
    REF_DBL area, max_area, normal[3], offset, center[3], hole[3];
    RSS(ref_export_faceid_range(ref_grid, &min_faceid, &max_faceid), "range");
    nhole = 0;
    for (faceid = min_faceid; faceid <= max_faceid; faceid++) {
      max_area = -1.0;
      largest_triangle = REF_EMPTY;
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        if (nodes[node_per] == faceid) {
          RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
          if (area > max_area) {
            max_area = area;
            largest_triangle = cell;
          }
        }
      }
      if (REF_EMPTY != largest_triangle) {
        RSS(ref_cell_nodes(ref_cell, largest_triangle, nodes), "tri nodes");
        RSS(ref_node_tri_normal(ref_node, nodes, normal), "normal");
        RSS(ref_math_normalize(normal), "norm");
        offset = 1.0e-4 * sqrt(max_area);
        if (offset > 1.0e-12) {
          nhole += 1;
        }
      }
    }
    fprintf(file, "%d\n", nhole);
    nhole = 0;
    for (faceid = min_faceid; faceid <= max_faceid; faceid++) {
      max_area = -1.0;
      largest_triangle = REF_EMPTY;
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        if (nodes[node_per] == faceid) {
          RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
          if (area > max_area) {
            max_area = area;
            largest_triangle = cell;
          }
        }
      }
      if (REF_EMPTY != largest_triangle) {
        RSS(ref_cell_nodes(ref_cell, largest_triangle, nodes), "tri nodes");
        RSS(ref_node_tri_normal(ref_node, nodes, normal), "normal");
        RSS(ref_math_normalize(normal), "norm");
        offset = 1.0e-4 * sqrt(max_area);
        if (offset > 1.0e-12) {
          RSS(ref_node_tri_centroid(ref_node, nodes, center), "center");
          hole[0] = center[0] - offset * normal[0];
          hole[1] = center[1] - offset * normal[1];
          hole[2] = center[2] - offset * normal[2];
          fprintf(file, "%d  %.16e  %.16e  %.16e\n", nhole, hole[0], hole[1],
                  hole[2]);
          nhole += 1;
        }
      }
    }
  }

  /* Part 4 - region attributes list */
  fprintf(file, "0\n"); /* no regions */

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_export_smesh(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nnode, ntri;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell;
  REF_INT dim, attr, mark;
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  nnode = ref_node_n(ref_node);

  dim = 3;
  attr = 0;
  mark = 0;

  fprintf(file, "%d  %d  %d  %d\n", nnode, dim, attr, mark);

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  for (node = 0; node < ref_node_n(ref_node); node++)
    fprintf(file, "%d  %.16e  %.16e  %.16e\n", node,
            ref_node_xyz(ref_node, 0, n2o[node]),
            ref_node_xyz(ref_node, 1, n2o[node]),
            ref_node_xyz(ref_node, 2, n2o[node]));

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  ntri = ref_cell_n(ref_cell);

  mark = 1;
  fprintf(file, "%d  %d\n", ntri, mark);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    fprintf(file, "%d  %d %d %d  %d\n", node_per, o2n[nodes[0]], o2n[nodes[1]],
            o2n[nodes[2]], nodes[3]);
  }

  fprintf(file, "0\n0\n");

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_export_fgrid(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nnode, ntri, ntet;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell;
  REF_INT ixyz;

  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  nnode = ref_node_n(ref_node);

  ntri = ref_cell_n(ref_grid_tri(ref_grid));

  ntet = ref_cell_n(ref_grid_tet(ref_grid));

  fprintf(file, "%d %d %d\n", nnode, ntri, ntet);

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  for (ixyz = 0; ixyz < 3; ixyz++)
    for (node = 0; node < ref_node_n(ref_node); node++)
      fprintf(file, " %.16e\n", ref_node_xyz(ref_node, ixyz, n2o[node]));

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < node_per; node++)
      fprintf(file, " %d", o2n[nodes[node]] + 1);
    fprintf(file, "\n");
  }
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    fprintf(file, " %d", nodes[3]);
    fprintf(file, "\n");
  }

  ref_cell = ref_grid_tet(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < node_per; node++)
      fprintf(file, " %d", o2n[nodes[node]] + 1);
    fprintf(file, "\n");
  }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

/* https://su2code.github.io/docs/Mesh-File/ */
static REF_STATUS ref_export_su2(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell;
  REF_INT group;
  REF_INT nnode;

  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  if (ref_grid_twod(ref_grid)) {
    fprintf(file, "NDIME= 2\n");
  } else {
    fprintf(file, "NDIME= 3\n");
  }

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  nnode = ref_node_n(ref_node);
  fprintf(file, "NPOIN= %d\n", nnode);

  if (ref_grid_twod(ref_grid)) {
    for (node = 0; node < ref_node_n(ref_node); node++)
      fprintf(file, " %.16e %.16e\n", ref_node_xyz(ref_node, 0, n2o[node]),
              ref_node_xyz(ref_node, 1, n2o[node]));
  } else {
    for (node = 0; node < ref_node_n(ref_node); node++)
      fprintf(file, " %.16e %.16e %.16e\n",
              ref_node_xyz(ref_node, 0, n2o[node]),
              ref_node_xyz(ref_node, 1, n2o[node]),
              ref_node_xyz(ref_node, 2, n2o[node]));
  }

  if (ref_grid_twod(ref_grid)) {
    REF_INT ntri = ref_cell_n(ref_grid_tri(ref_grid));
    REF_INT nqua = ref_cell_n(ref_grid_qua(ref_grid));

    fprintf(file, "NELEM= %d\n", ntri + nqua);
    /* uses VTK windings */
    each_ref_grid_2d_ref_cell(ref_grid, group, ref_cell) {
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        switch (ref_cell_type(ref_cell)) {
          case REF_CELL_TRI:
            fprintf(file, "5");
            break;
          case REF_CELL_QUA:
            fprintf(file, "9");
            break;
          case REF_CELL_EDG:
          case REF_CELL_ED2:
          case REF_CELL_ED3:
          case REF_CELL_TR2:
          case REF_CELL_TR3:
          case REF_CELL_TET:
          case REF_CELL_PYR:
          case REF_CELL_PRI:
          case REF_CELL_HEX:
            RSS(REF_IMPLEMENT, "2D SU2 element");
            break;
        }
        for (node = 0; node < ref_cell_node_per(ref_cell); node++)
          fprintf(file, " %d", o2n[nodes[node]]);
        fprintf(file, "\n");
      }
    }
  } else {
    REF_INT ntet = ref_cell_n(ref_grid_tet(ref_grid));
    REF_INT npyr = ref_cell_n(ref_grid_pyr(ref_grid));
    REF_INT npri = ref_cell_n(ref_grid_pri(ref_grid));
    REF_INT nhex = ref_cell_n(ref_grid_hex(ref_grid));

    fprintf(file, "NELEM= %d\n", ntet + npyr + npri + nhex);
    /* uses VTK windings */
    each_ref_grid_3d_ref_cell(ref_grid, group, ref_cell) {
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        switch (ref_cell_type(ref_cell)) {
          case REF_CELL_TET:
            fprintf(file, "10");
            break;
          case REF_CELL_PYR:
            fprintf(file, "14");
            VTK_PYRAMID_ORDER(nodes);
            break;
          case REF_CELL_PRI:
            fprintf(file, "13");
            VTK_WEDGE_ORDER(nodes);
            break;
          case REF_CELL_HEX:
            fprintf(file, "12");
            break;
          case REF_CELL_EDG:
          case REF_CELL_ED2:
          case REF_CELL_ED3:
          case REF_CELL_TRI:
          case REF_CELL_TR2:
          case REF_CELL_TR3:
          case REF_CELL_QUA:
            RSS(REF_IMPLEMENT, "3D SU2 element");
            break;
        }
        for (node = 0; node < ref_cell_node_per(ref_cell); node++)
          fprintf(file, " %d", o2n[nodes[node]]);
        fprintf(file, "\n");
      }
    }
  }

  if (ref_grid_twod(ref_grid)) {
    REF_INT faceid, min_faceid, max_faceid;
    REF_INT nedg;
    RSS(ref_export_cell_id_range(ref_grid_edg(ref_grid), &min_faceid,
                                 &max_faceid),
        "range");
    fprintf(file, "NMARK= %d\n", max_faceid - min_faceid + 1);

    for (faceid = min_faceid; faceid <= max_faceid; faceid++) {
      ref_cell = ref_grid_edg(ref_grid);
      node_per = ref_cell_node_per(ref_cell);
      nedg = 0;
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        if (nodes[node_per] == faceid) {
          nedg++;
        }
      }
      fprintf(file, "MARKER_TAG= %d\n", faceid);
      fprintf(file, "MARKER_ELEMS= %d\n", nedg);

      ref_cell = ref_grid_edg(ref_grid);
      node_per = ref_cell_node_per(ref_cell);
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        if (nodes[node_per] == faceid) {
          fprintf(file, "3");
          for (node = 0; node < node_per; node++)
            fprintf(file, " %d", o2n[nodes[node]]);
          fprintf(file, "\n");
        }
      }
    }

  } else {
    REF_INT faceid, min_faceid, max_faceid;
    REF_INT ntri;
    REF_INT nqua;
    RSS(ref_export_faceid_range(ref_grid, &min_faceid, &max_faceid), "range");
    fprintf(file, "NMARK= %d\n", max_faceid - min_faceid + 1);

    for (faceid = min_faceid; faceid <= max_faceid; faceid++) {
      ref_cell = ref_grid_tri(ref_grid);
      node_per = ref_cell_node_per(ref_cell);
      ntri = 0;
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        if (nodes[node_per] == faceid) {
          ntri++;
        }
      }

      ref_cell = ref_grid_qua(ref_grid);
      node_per = ref_cell_node_per(ref_cell);
      nqua = 0;
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        if (nodes[node_per] == faceid) {
          nqua++;
        }
      }

      fprintf(file, "MARKER_TAG= %d\n", faceid);
      fprintf(file, "MARKER_ELEMS= %d\n", ntri + nqua);

      ref_cell = ref_grid_tri(ref_grid);
      node_per = ref_cell_node_per(ref_cell);
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        if (nodes[node_per] == faceid) {
          fprintf(file, "5");
          for (node = 0; node < node_per; node++)
            fprintf(file, " %d", o2n[nodes[node]]);
          fprintf(file, "\n");
        }
      }

      ref_cell = ref_grid_qua(ref_grid);
      node_per = ref_cell_node_per(ref_cell);
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        if (nodes[node_per] == faceid) {
          fprintf(file, "9");
          for (node = 0; node < node_per; node++)
            fprintf(file, " %d", o2n[nodes[node]]);
          fprintf(file, "\n");
        }
      }
    }
  }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_export_ugrid(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell;
  REF_INT group;
  REF_INT faceid, min_faceid, max_faceid;

  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  nnode = ref_node_n(ref_node);

  ntri = ref_cell_n(ref_grid_tri(ref_grid));
  nqua = ref_cell_n(ref_grid_qua(ref_grid));

  ntet = ref_cell_n(ref_grid_tet(ref_grid));
  npyr = ref_cell_n(ref_grid_pyr(ref_grid));
  npri = ref_cell_n(ref_grid_pri(ref_grid));
  nhex = ref_cell_n(ref_grid_hex(ref_grid));

  fprintf(file, "%d %d %d %d %d %d %d\n", nnode, ntri, nqua, ntet, npyr, npri,
          nhex);

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  for (node = 0; node < ref_node_n(ref_node); node++)
    fprintf(file, " %.16e %.16e %.16e\n", ref_node_xyz(ref_node, 0, n2o[node]),
            ref_node_xyz(ref_node, 1, n2o[node]),
            ref_node_xyz(ref_node, 2, n2o[node]));

  RSS(ref_export_faceid_range(ref_grid, &min_faceid, &max_faceid), "range");

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for (faceid = min_faceid; faceid <= max_faceid; faceid++)
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[node_per] == faceid) {
        for (node = 0; node < node_per; node++)
          fprintf(file, " %d", o2n[nodes[node]] + 1);
        fprintf(file, "\n");
      }
    }

  ref_cell = ref_grid_qua(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for (faceid = min_faceid; faceid <= max_faceid; faceid++)
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[node_per] == faceid) {
        for (node = 0; node < node_per; node++)
          fprintf(file, " %d", o2n[nodes[node]] + 1);
        fprintf(file, "\n");
      }
    }

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for (faceid = min_faceid; faceid <= max_faceid; faceid++) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[node_per] == faceid) {
        fprintf(file, " %d", nodes[3]);
        fprintf(file, "\n");
      }
    }
  }

  ref_cell = ref_grid_qua(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for (faceid = min_faceid; faceid <= max_faceid; faceid++) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[node_per] == faceid) {
        fprintf(file, " %d", nodes[4]);
        fprintf(file, "\n");
      }
    }
  }

  each_ref_grid_3d_ref_cell(ref_grid, group, ref_cell) {
    node_per = ref_cell_node_per(ref_cell);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      for (node = 0; node < node_per; node++)
        fprintf(file, " %d", o2n[nodes[node]] + 1);
      fprintf(file, "\n");
    }
  }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_export_bin_ugrid_int(FILE *file, REF_BOOL swap,
                                           REF_BOOL fat, REF_INT output) {
  if (fat) {
    REF_LONG actual;
    actual = (REF_LONG)output;
    if (swap) SWAP_LONG(actual);
    REIS(1, fwrite(&actual, sizeof(REF_LONG), 1, file), "output long");
  } else {
    REF_INT actual;
    actual = output;
    if (swap) SWAP_INT(actual);
    REIS(1, fwrite(&actual, sizeof(REF_INT), 1, file), "output int");
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_export_bin_ugrid(REF_GRID ref_grid, const char *filename,
                                       REF_BOOL swap, REF_BOOL fat) {
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per, cell;
  REF_DBL swapped_dbl;
  REF_INT group;
  REF_INT faceid, min_faceid, max_faceid;

  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RSS(ref_export_bin_ugrid_int(file, swap, fat, ref_node_n(ref_node)), "nnode");
  RSS(ref_export_bin_ugrid_int(file, swap, fat,
                               ref_cell_n(ref_grid_tri(ref_grid))),
      "ntri");
  RSS(ref_export_bin_ugrid_int(file, swap, fat,
                               ref_cell_n(ref_grid_qua(ref_grid))),
      "nqua");
  RSS(ref_export_bin_ugrid_int(file, swap, fat,
                               ref_cell_n(ref_grid_tet(ref_grid))),
      "ntet");
  RSS(ref_export_bin_ugrid_int(file, swap, fat,
                               ref_cell_n(ref_grid_pyr(ref_grid))),
      "npyr");
  RSS(ref_export_bin_ugrid_int(file, swap, fat,
                               ref_cell_n(ref_grid_pri(ref_grid))),
      "npri");
  RSS(ref_export_bin_ugrid_int(file, swap, fat,
                               ref_cell_n(ref_grid_hex(ref_grid))),
      "nhex");

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  for (node = 0; node < ref_node_n(ref_node); node++) {
    swapped_dbl = ref_node_xyz(ref_node, 0, n2o[node]);
    if (swap) SWAP_DBL(swapped_dbl);
    REIS(1, fwrite(&swapped_dbl, sizeof(REF_DBL), 1, file), "x");
    swapped_dbl = ref_node_xyz(ref_node, 1, n2o[node]);
    if (swap) SWAP_DBL(swapped_dbl);
    REIS(1, fwrite(&swapped_dbl, sizeof(REF_DBL), 1, file), "y");
    swapped_dbl = ref_node_xyz(ref_node, 2, n2o[node]);
    if (swap) SWAP_DBL(swapped_dbl);
    REIS(1, fwrite(&swapped_dbl, sizeof(REF_DBL), 1, file), "z");
  }

  RSS(ref_export_faceid_range(ref_grid, &min_faceid, &max_faceid), "range");

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for (faceid = min_faceid; faceid <= max_faceid; faceid++) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[node_per] == faceid) {
        for (node = 0; node < node_per; node++) {
          nodes[node] = o2n[nodes[node]] + 1;
          RSS(ref_export_bin_ugrid_int(file, swap, fat, nodes[node]), "c2n");
        }
      }
    }
  }

  ref_cell = ref_grid_qua(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for (faceid = min_faceid; faceid <= max_faceid; faceid++) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[node_per] == faceid) {
        for (node = 0; node < node_per; node++) {
          nodes[node] = o2n[nodes[node]] + 1;
          RSS(ref_export_bin_ugrid_int(file, swap, fat, nodes[node]), "c2n");
        }
      }
    }
  }

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for (faceid = min_faceid; faceid <= max_faceid; faceid++) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[node_per] == faceid) {
        RSS(ref_export_bin_ugrid_int(file, swap, fat, nodes[3]), "c2n");
      }
    }
  }

  ref_cell = ref_grid_qua(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for (faceid = min_faceid; faceid <= max_faceid; faceid++) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[node_per] == faceid) {
        RSS(ref_export_bin_ugrid_int(file, swap, fat, nodes[4]), "c2n");
      }
    }
  }

  each_ref_grid_3d_ref_cell(ref_grid, group, ref_cell) {
    node_per = ref_cell_node_per(ref_cell);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      for (node = 0; node < node_per; node++) {
        nodes[node] = o2n[nodes[node]] + 1;
        RSS(ref_export_bin_ugrid_int(file, swap, fat, nodes[node]), "c2n");
      }
    }
  }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_export_c(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT node, ixyz;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;

  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  for (node = 0; node < ref_node_n(ref_node); node++) {
    fprintf(file, "  RSS(ref_node_add(ref_node,%d,&node),\"node\");\n", node);
    for (ixyz = 0; ixyz < 3; ixyz++)
      fprintf(file, "  ref_node_xyz(ref_node,%d,node) = %.15e;\n", ixyz,
              ref_node_xyz(ref_node, ixyz, n2o[node]));
  }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < ref_cell_size_per(ref_cell); node++)
      fprintf(file, "  nodes[%d] = %d;\n", node, o2n[nodes[node]]);
    fprintf(
        file,
        "  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),\"qua\");\n");
  }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_export_eps(REF_GRID ref_grid, const char *filename) {
  FILE *f;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT node, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  f = fopen("ref_export.gnuplot", "w");
  if (NULL == (void *)f) printf("unable to open ref_export.gnuplot\n");
  RNS(f, "unable to open file");

  fprintf(f, "reset\n");
  fprintf(f, "set term postscript eps\n");
  fprintf(f, "set output '%s'\n", filename);
  fprintf(f, "set size ratio -1\n");
  fprintf(f, "set xlabel 'X'\n");
  fprintf(f, "set ylabel 'Z'\n");
  fprintf(f, "plot '-' title '' with lines lw 0.5\n");

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      fprintf(f, "%25.15f %25.15f\n", ref_node_xyz(ref_node, 0, nodes[node]),
              ref_node_xyz(ref_node, 2, nodes[node]));
    }
    node = 0;
    fprintf(f, "%25.15f %25.15f\n", ref_node_xyz(ref_node, 0, nodes[node]),
            ref_node_xyz(ref_node, 2, nodes[node]));
    fprintf(f, "\n\n");
  }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      fprintf(f, "%25.15f %25.15f\n", ref_node_xyz(ref_node, 0, nodes[node]),
              ref_node_xyz(ref_node, 2, nodes[node]));
    }
    node = 0;
    fprintf(f, "%25.15f %25.15f\n", ref_node_xyz(ref_node, 0, nodes[node]),
            ref_node_xyz(ref_node, 2, nodes[node]));
    fprintf(f, "\n\n");
  }

  fprintf(f, "e\n");
  fclose(f);

  REIS(0, system("gnuplot ref_export.gnuplot"), "gnuplot failed");
  REIS(0, remove("ref_export.gnuplot"), "temp clean up");

  return REF_SUCCESS;
}

static REF_STATUS ref_export_pdf(REF_GRID ref_grid, const char *pdf_filename) {
  char temp_filename[] = "ref_export_temp_for_pdf.eps";
  char command[1024];
  RSS(ref_export_eps(ref_grid, temp_filename), "temp eps");
  sprintf(command, "epstopdf %s -o=%s", temp_filename, pdf_filename);
  REIS(0, system(command), "epstopdf failed");
  REIS(0, remove(temp_filename), "temp clean up");

  return REF_SUCCESS;
}

static REF_STATUS ref_export_html(REF_GRID ref_grid, const char *filename) {
  FILE *f;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT node, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT *o2n, *n2o;

  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  fprintf(f, "<html>\n");
  fprintf(f, "  <head>\n");
  fprintf(f, "    <title>refine export</title>\n");
  fprintf(f, "    <link rel='stylesheet' type='text/css'\n");
  fprintf(f, "          href='http://www.x3dom.org/download/x3dom.css'>\n");
  fprintf(f, "    </link>\n");
  fprintf(f, "    <script type='text/javascript'\n");
  fprintf(f, "          src='http://www.x3dom.org/download/x3dom.js'>\n");
  fprintf(f, "    </script>\n");
  fprintf(f, "    <style>\n");
  fprintf(f, "      x3d {width:100%%;height:100%%;border:none}\n");
  fprintf(f, "      body {margin:0;width:100%%;height:100%%;}\n");
  fprintf(f, "    </style>\n");
  fprintf(f, "  </head>\n");
  fprintf(f, "  <body id='body'>\n");
  fprintf(f, "    <a href=\"http://x3dom.org/docs/dev/navigation.html\">\n");
  fprintf(f, "camera control help\n");
  fprintf(f, "    </a>\n");
  fprintf(f, "    <x3d id='x3d'><scene><shape>\n");

  fprintf(f, "      <IndexedLineSet coordIndex='\n");
  if (REF_TRUE) {
    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      for (node = 0; node < ref_cell_node_per(ref_cell); node++)
        fprintf(f, " %d", o2n[nodes[node]]);
      fprintf(f, " %d %d\n", o2n[nodes[0]], -1);
    }
  }

  if (REF_TRUE) {
    ref_cell = ref_grid_qua(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      for (node = 0; node < ref_cell_node_per(ref_cell); node++)
        fprintf(f, " %d", o2n[nodes[node]]);
      fprintf(f, " %d %d\n", o2n[nodes[0]], -1);
    }
  }

  if (REF_TRUE) {
    ref_cell = ref_grid_tet(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      fprintf(f, " %d %d %d %d %d %d %d %d\n", o2n[nodes[0]], o2n[nodes[1]],
              o2n[nodes[2]], o2n[nodes[3]], o2n[nodes[0]], o2n[nodes[3]],
              o2n[nodes[1]], -1);
    }
  }

  fprintf(f, "      ' >\n");

  fprintf(f, "      <Coordinate point='\n");
  for (node = 0; node < ref_node_n(ref_node); node++) {
    fprintf(f, " %.15e", ref_node_xyz(ref_node, 0, n2o[node]));
    fprintf(f, " %.15e", ref_node_xyz(ref_node, 1, n2o[node]));
    fprintf(f, " %.15e\n", ref_node_xyz(ref_node, 2, n2o[node]));
  }
  fprintf(f, "      ' />\n");

  fprintf(f, "      </IndexedLineSet>\n");

  fprintf(f, "    </shape></scene></x3d>\n");
  fprintf(f, "  </body>\n");
  fprintf(f, "</html>\n");

  ref_free(n2o);
  ref_free(o2n);

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_export_meshb_next_position(FILE *file, REF_INT version,
                                          REF_FILEPOS next_position) {
  int32_t one_word;
  int64_t two_word;

  if (3 <= version) {
    two_word = (int64_t)next_position;
    REIS(1, fwrite(&two_word, sizeof(two_word), 1, file), "write next pos");
  } else {
    if (next_position < -2147483647 || 2147483647 < next_position) {
      printf("next_position %ld outside int32 limits %d %d\n",
             (REF_LONG)next_position, -2147483647, 2147483647);
      RSS(REF_INVALID, "meshb version does not support file size");
    }
    one_word = (int32_t)next_position;
    REIS(1, fwrite(&one_word, sizeof(one_word), 1, file), "write next pos");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_export_meshb_int(FILE *file, REF_INT version,
                                       REF_INT value) {
  int int_value;
  long long_value;
  if (version < 4) {
    int_value = (int)value;
    REIS(1, fwrite(&int_value, sizeof(int), 1, file), "int value");
  } else {
    long_value = (long)value;
    REIS(1, fwrite(&long_value, sizeof(long), 1, file), "long value");
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_export_meshb(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell;
  REF_INT *o2n, *n2o;
  REF_INT code, version, dim;
  REF_FILEPOS next_position;
  REF_INT keyword_code, header_size, int_size, fp_size;
  REF_INT node;
  REF_INT group, node_per, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER + 1]; /* everyone gets id in meshb */
  REF_INT type, id, i, size_bytes;
  REF_INT n0, n1, n2, n3, n4;
  REF_INT geom;
  int ngeom;

  dim = 3;
  if (ref_grid_twod(ref_grid)) dim = 2;

  version = 2;
  if (1 < ref_grid_meshb_version(ref_grid)) {
    version = ref_grid_meshb_version(ref_grid);
  } else {
    if (REF_EXPORT_MESHB_VERTEX_3 < ref_node_n_global(ref_node)) version = 3;
    if (REF_EXPORT_MESHB_VERTEX_4 < ref_node_n_global(ref_node)) version = 4;
  }

  int_size = 4;
  fp_size = 4;
  if (2 < version) fp_size = 8;
  if (3 < version) int_size = 8;
  header_size = 4 + fp_size + int_size;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  code = 1;
  REIS(1, fwrite(&code, sizeof(int), 1, file), "code");
  REIS(1, fwrite(&version, sizeof(int), 1, file), "version");
  next_position = (REF_FILEPOS)(4 + fp_size + 4) + ftell(file);
  keyword_code = 3;
  REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "dim code");
  RSS(ref_export_meshb_next_position(file, version, next_position), "next pos");
  REIS(1, fwrite(&dim, sizeof(int), 1, file), "dim");
  REIS(next_position, ftell(file), "dim inconsistent");

  if (ref_node_n(ref_node) > 0) {
    next_position =
        (REF_FILEPOS)header_size +
        (REF_FILEPOS)ref_node_n(ref_node) * (REF_FILEPOS)(dim * 8 + int_size) +
        ftell(file);
    keyword_code = 4;
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "vertex version code");
    RSS(ref_export_meshb_next_position(file, version, next_position), "next p");
    RSS(ref_export_meshb_int(file, version, ref_node_n(ref_node)), "nnode");
    for (node = 0; node < ref_node_n(ref_node); node++) {
      REIS(1,
           fwrite(&(ref_node_xyz(ref_node, 0, n2o[node])), sizeof(double), 1,
                  file),
           "x");
      REIS(1,
           fwrite(&(ref_node_xyz(ref_node, 1, n2o[node])), sizeof(double), 1,
                  file),
           "y");
      if (!ref_grid_twod(ref_grid))
        REIS(1,
             fwrite(&(ref_node_xyz(ref_node, 2, n2o[node])), sizeof(double), 1,
                    file),
             "z");
      RSS(ref_export_meshb_int(file, version, REF_EXPORT_MESHB_VERTEX_ID),
          "nnode");
    }
    REIS(next_position, ftell(file), "vertex inconsistent");
  }

  each_ref_grid_all_ref_cell(ref_grid, group, ref_cell) {
    if (ref_cell_n(ref_cell) > 0) {
      RSS(ref_cell_meshb_keyword(ref_cell, &keyword_code), "kw");
      node_per = ref_cell_node_per(ref_cell);
      next_position = ftell(file) + (REF_FILEPOS)header_size +
                      (REF_FILEPOS)ref_cell_n(ref_cell) *
                          (REF_FILEPOS)(int_size * (node_per + 1));
      REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "keyword code");
      RSS(ref_export_meshb_next_position(file, version, next_position), "next");
      RSS(ref_export_meshb_int(file, version, ref_cell_n(ref_cell)), "ncell");
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        for (node = 0; node < node_per; node++) {
          nodes[node] = o2n[nodes[node]] + 1;
        }
        if (REF_CELL_PYR == ref_cell_type(ref_cell)) {
          /* convention: square basis is 0-1-2-3
             (oriented counter clockwise like trias) and top vertex is 4 */
          n0 = nodes[0];
          n1 = nodes[3];
          n2 = nodes[4];
          n3 = nodes[1];
          n4 = nodes[2];
          nodes[0] = n0;
          nodes[1] = n1;
          nodes[2] = n2;
          nodes[3] = n3;
          nodes[4] = n4;
        }
        if (!ref_cell_last_node_is_an_id(ref_cell))
          nodes[node_per] = REF_EXPORT_MESHB_3D_ID;
        for (node = 0; node < (1 + node_per); node++) {
          RSS(ref_export_meshb_int(file, version, nodes[node]), "c2n");
        }
      }
      REIS(next_position, ftell(file), "cell inconsistent");
    }
  }

  each_ref_type(ref_geom, type) {
    ngeom = 0;
    each_ref_geom_of(ref_geom, type, geom) ngeom++;
    if (ngeom > 0) {
      keyword_code = 40 + type; /* GmfVerticesOnGeometricVertices */
      next_position =
          ftell(file) + (REF_FILEPOS)header_size +
          (REF_FILEPOS)ngeom *
              (REF_FILEPOS)(int_size * 2 + 8 * type + (0 < type ? 8 : 0));
      REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "keyword");
      RSS(ref_export_meshb_next_position(file, version, next_position), "next");
      RSS(ref_export_meshb_int(file, version, ngeom), "ngeom");
      each_ref_geom_of(ref_geom, type, geom) {
        node = o2n[ref_geom_node(ref_geom, geom)] + 1;
        id = ref_geom_id(ref_geom, geom);
        RSS(ref_export_meshb_int(file, version, node), "node");
        RSS(ref_export_meshb_int(file, version, id), "id");
        for (i = 0; i < type; i++)
          REIS(1,
               fwrite(&(ref_geom_param(ref_geom, i, geom)), sizeof(double), 1,
                      file),
               "id");
        if (0 < type) {
          double double_gref = (double)ref_geom_gref(ref_geom, geom);
          REIS(1, fwrite(&(double_gref), sizeof(double), 1, file), "gref");
        }
      }
      REIS(next_position, ftell(file), "geom inconsistent");
    }
  }

  if (0 < ref_geom_cad_data_size(ref_geom)) {
    keyword_code = 126; /* GmfByteFlow 173-47 */
    next_position = (REF_FILEPOS)header_size +
                    (REF_FILEPOS)ref_geom_cad_data_size(ref_geom) + ftell(file);
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "keyword");
    RSS(ref_export_meshb_next_position(file, version, next_position), "next p");
    size_bytes = (REF_INT)ref_geom_cad_data_size(ref_geom);
    RSS(ref_export_meshb_int(file, version, size_bytes), "size in bytes");
    REIS(ref_geom_cad_data_size(ref_geom),
         fwrite(ref_geom_cad_data(ref_geom), sizeof(REF_BYTE),
                ref_geom_cad_data_size(ref_geom), file),
         "node");
    REIS(next_position, ftell(file), "cad_model inconsistent");
  }

  /* End */
  keyword_code = 54; /* GmfEnd 101-47 */
  REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "vertex version code");
  next_position = 0;
  RSS(ref_export_meshb_next_position(file, version, next_position), "next p");

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_export_twod_msh(REF_GRID ref_grid, const char *filename) {
  FILE *f;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nnode;
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT nedge;
  REF_INT ntri;

  RAS(ref_grid_twod(ref_grid), "expected twod convention grid");

  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  fprintf(f, "MeshVersionFormatted 0\n\n");
  fprintf(f, "Dimension 2\n\n");

  ref_malloc_init(o2n, ref_node_max(ref_node), REF_INT, REF_EMPTY);
  ref_malloc_init(n2o, ref_node_max(ref_node), REF_INT, REF_EMPTY);

  nnode = 0;
  each_ref_node_valid_node(ref_node, node) {
    o2n[node] = nnode;
    n2o[nnode] = node;
    nnode++;
  }

  fprintf(f, "\nVertices\n%d\n", nnode);
  for (node = 0; node < nnode; node++) {
    fprintf(f, "%.16E %.16E %d\n", ref_node_xyz(ref_node, 0, n2o[node]),
            ref_node_xyz(ref_node, 2, n2o[node]), 1);
  }

  ref_cell = ref_grid_edg(ref_grid);
  fprintf(f, "\nEdges\n%d\n", ref_cell_n(ref_cell));
  nedge = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    nedge++;
    fprintf(f, "%d %d %d\n", o2n[nodes[0]] + 1, o2n[nodes[1]] + 1, nodes[2]);
  }
  REIS(nedge, ref_cell_n(ref_cell), "edge/quad miscount");

  ref_cell = ref_grid_tri(ref_grid);
  fprintf(f, "\nTriangles\n%d\n", ref_cell_n(ref_cell));
  ntri = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    ntri++;
    fprintf(f, "%d %d %d %d\n", o2n[nodes[0]] + 1, o2n[nodes[2]] + 1,
            o2n[nodes[1]] + 1, nodes[3]);
  }
  REIS(ntri, ref_cell_n(ref_cell), "triangle miscount");

  ref_free(n2o);
  ref_free(o2n);

  fclose(f);

  return REF_SUCCESS;
}

static REF_STATUS ref_export_i_like_cfd_grid(REF_GRID ref_grid,
                                             const char *filename) {
  FILE *f;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_INT *o2n, *n2o, *c2n, *order;
  REF_INT nnode;
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT min_id, max_id, id;
  REF_INT ntri, nquad, nedge, edge;

  RAS(ref_grid_twod(ref_grid), "expected twod convention grid");

  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  ref_malloc_init(o2n, ref_node_max(ref_node), REF_INT, REF_EMPTY);
  ref_malloc_init(n2o, ref_node_max(ref_node), REF_INT, REF_EMPTY);

  nnode = 0;
  each_ref_node_valid_node(ref_node, node) {
    o2n[node] = nnode;
    n2o[nnode] = node;
    nnode++;
  }

  fprintf(f, "%d %d %d\n", nnode, ref_cell_n(ref_grid_tri(ref_grid)),
          ref_cell_n(ref_grid_qua(ref_grid)));

  for (node = 0; node < nnode; node++) {
    fprintf(f, "%.16E %.16E\n", ref_node_xyz(ref_node, 0, n2o[node]),
            ref_node_xyz(ref_node, 1, n2o[node]));
  }

  ref_cell = ref_grid_tri(ref_grid);
  ntri = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    ntri++;
    fprintf(f, "%d %d %d\n", o2n[nodes[0]] + 1, o2n[nodes[1]] + 1,
            o2n[nodes[2]] + 1);
  }
  REIS(ntri, ref_cell_n(ref_cell), "triangle miscount");

  ref_cell = ref_grid_qua(ref_grid);
  nquad = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    ntri++;
    fprintf(f, "%d %d %d %d\n", o2n[nodes[0]] + 1, o2n[nodes[1]] + 1,
            o2n[nodes[2]] + 1, o2n[nodes[3]] + 1);
  }
  REIS(nquad, ref_cell_n(ref_cell), "quad miscount");

  ref_cell = ref_grid_edg(ref_grid);
  RSS(ref_cell_id_range(ref_cell, ref_grid_mpi(ref_grid), &min_id, &max_id),
      "id range");
  fprintf(f, "%d\n", max_id - min_id + 1);
  for (id = min_id; id <= max_id; id++) {
    nedge = 0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[2] == id) {
        nedge++;
      }
    }
    fprintf(f, "%d\n", nedge + 1); /* +1, both end points */
  }

  ref_malloc(c2n, 2 * ref_cell_n(ref_cell), REF_INT);
  ref_malloc(order, ref_cell_n(ref_cell), REF_INT);
  for (id = min_id; id <= max_id; id++) {
    nedge = 0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[2] == id) {
        RSS(ref_grid_orient_edg(ref_grid, nodes), "orient based on tri");
        c2n[0 + 2 * nedge] = nodes[0];
        c2n[1 + 2 * nedge] = nodes[1];
        nedge++;
      }
    }
    if (nedge > 0) {
      RSS(ref_export_order_segments(nedge, c2n, order), "order");
      fprintf(f, "%d\n", o2n[c2n[0 + 2 * order[0]]] + 1);
      for (edge = 0; edge < nedge; edge++) {
        fprintf(f, "%d\n", o2n[c2n[1 + 2 * order[edge]]] + 1);
      }
    }
  }

  ref_free(order);
  ref_free(c2n);
  ref_free(n2o);
  ref_free(o2n);

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_export_by_extension(REF_GRID ref_grid, const char *filename) {
  size_t end_of_string;

  end_of_string = strlen(filename);

  if (strcmp(&filename[end_of_string - 4], ".vtk") == 0) {
    RSS(ref_export_vtk(ref_grid, filename), "vtk export failed");
  } else if (strcmp(&filename[end_of_string - 2], ".c") == 0) {
    RSS(ref_export_c(ref_grid, filename), "C export failed");
  } else if (strcmp(&filename[end_of_string - 4], ".tec") == 0) {
    RSS(ref_export_tec(ref_grid, filename), "tec export failed");
  } else if (strcmp(&filename[end_of_string - 4], ".plt") == 0) {
    RSS(ref_gather_by_extension(ref_grid, filename), "plt gather failed");
  } else if (strcmp(&filename[end_of_string - 4], ".eps") == 0) {
    RSS(ref_export_eps(ref_grid, filename), "eps export failed");
  } else if (strcmp(&filename[end_of_string - 4], ".pdf") == 0) {
    RSS(ref_export_pdf(ref_grid, filename), "pdf export failed");
  } else if (strcmp(&filename[end_of_string - 4], ".su2") == 0) {
    RSS(ref_export_su2(ref_grid, filename), "su2 export failed");
  } else if (strcmp(&filename[end_of_string - 10], ".lb8.ugrid") == 0) {
    RSS(ref_export_bin_ugrid(ref_grid, filename, REF_FALSE, REF_FALSE),
        "lb8.ugrid export failed");
  } else if (strcmp(&filename[end_of_string - 9], ".b8.ugrid") == 0) {
    RSS(ref_export_bin_ugrid(ref_grid, filename, REF_TRUE, REF_FALSE),
        "b8.ugrid export failed");
  } else if (strcmp(&filename[end_of_string - 11], ".lb8l.ugrid") == 0) {
    RSS(ref_export_bin_ugrid(ref_grid, filename, REF_FALSE, REF_TRUE),
        "lb8l.ugrid export failed");
  } else if (strcmp(&filename[end_of_string - 10], ".b8l.ugrid") == 0) {
    RSS(ref_export_bin_ugrid(ref_grid, filename, REF_TRUE, REF_TRUE),
        "b8l.ugrid export failed");
  } else if (strcmp(&filename[end_of_string - 12], ".lb8.ugrid64") == 0) {
    RSS(ref_export_bin_ugrid(ref_grid, filename, REF_FALSE, REF_TRUE),
        "lb8.ugrid64 export failed");
  } else if (strcmp(&filename[end_of_string - 11], ".b8.ugrid64") == 0) {
    RSS(ref_export_bin_ugrid(ref_grid, filename, REF_TRUE, REF_TRUE),
        "b8.ugrid64 export failed");
  } else if (strcmp(&filename[end_of_string - 6], ".ugrid") == 0) {
    RSS(ref_export_ugrid(ref_grid, filename), "ugrid export failed");
  } else if (strcmp(&filename[end_of_string - 5], ".grid") == 0) {
    RSS(ref_export_i_like_cfd_grid(ref_grid, filename), "grid export failed");
  } else if (strcmp(&filename[end_of_string - 5], ".poly") == 0) {
    RSS(ref_export_poly(ref_grid, filename), "poly export failed");
  } else if (strcmp(&filename[end_of_string - 6], ".smesh") == 0) {
    RSS(ref_export_smesh(ref_grid, filename), "smesh export failed");
  } else if (strcmp(&filename[end_of_string - 6], ".fgrid") == 0) {
    RSS(ref_export_fgrid(ref_grid, filename), "fgrid export failed");
  } else if (strcmp(&filename[end_of_string - 5], ".html") == 0) {
    RSS(ref_export_html(ref_grid, filename), "html export failed");
  } else if (strcmp(&filename[end_of_string - 6], ".meshb") == 0) {
    RSS(ref_export_meshb(ref_grid, filename), "meshb export failed");
  } else if (strcmp(&filename[end_of_string - 4], ".msh") == 0) {
    RSS(ref_export_twod_msh(ref_grid, filename), "msh export failed");
  } else {
    RSS(ref_gather_by_extension(ref_grid, filename), "export via gather");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_export_order_segments(REF_INT n, REF_INT *c2n, REF_INT *order) {
  REF_INT i, j, first_cell;
  REF_BOOL found, loop;

  if (0 == n) return REF_SUCCESS;

  first_cell = REF_EMPTY;
  for (i = 0; i < n; i++) {
    found = REF_FALSE;
    for (j = 0; j < n; j++) {
      if (c2n[0 + 2 * i] == c2n[1 + 2 * j]) {
        found = REF_TRUE;
        break;
      }
    }
    if (!found) {
      first_cell = i;
      break;
    }
  }

  /* arbitrary for loops */
  loop = REF_FALSE;
  if (REF_EMPTY == first_cell) {
    loop = REF_TRUE;
    first_cell = 0;
  }

  order[0] = first_cell;
  for (i = 1; i < n; i++) {
    found = REF_FALSE;
    for (j = 0; j < n; j++) {
      if (c2n[1 + 2 * order[i - 1]] == c2n[0 + 2 * j]) {
        order[i] = j;
        found = REF_TRUE;
        break;
      }
    }
    RAS(found, "next segment not found");
  }

  if (loop)
    REIS(c2n[0 + 2 * order[0]], c2n[1 + 2 * order[n - 1]],
         "first and last do not match for loop");

  return REF_SUCCESS;
}
