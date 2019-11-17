
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

/*
  tecplot "brick"
      7---6
     /|  /|
    4-+-5 |
    | | | |
    | 3-+-2
    |/  |/
    0---1
 */

#define TEC_BRICK_TET(brick, nodes) \
  {                                 \
    brick[0] = nodes[0];            \
    brick[1] = nodes[1];            \
    brick[2] = nodes[2];            \
    brick[3] = nodes[2];            \
    brick[4] = nodes[3];            \
    brick[5] = nodes[3];            \
    brick[6] = nodes[3];            \
    brick[7] = nodes[3];            \
  }

#define TEC_BRICK_PYR(brick, nodes) \
  {                                 \
    brick[0] = nodes[0];            \
    brick[1] = nodes[1];            \
    brick[2] = nodes[2];            \
    brick[3] = nodes[2];            \
    brick[4] = nodes[3];            \
    brick[5] = nodes[4];            \
    brick[6] = nodes[2];            \
    brick[7] = nodes[2];            \
  }

#define TEC_BRICK_PRI(brick, nodes) \
  {                                 \
    brick[0] = nodes[0];            \
    brick[1] = nodes[1];            \
    brick[2] = nodes[2];            \
    brick[3] = nodes[2];            \
    brick[4] = nodes[3];            \
    brick[5] = nodes[4];            \
    brick[6] = nodes[5];            \
    brick[7] = nodes[5];            \
  }

#define TEC_BRICK_HEX(brick, nodes) \
  {                                 \
    brick[0] = nodes[0];            \
    brick[1] = nodes[1];            \
    brick[2] = nodes[2];            \
    brick[3] = nodes[3];            \
    brick[4] = nodes[4];            \
    brick[5] = nodes[5];            \
    brick[6] = nodes[6];            \
    brick[7] = nodes[7];            \
  }

/* https://www.vtk.org/VTK/img/file-formats.pdf */
REF_STATUS ref_export_vtk(REF_GRID ref_grid, const char *filename) {
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

  ref_cell = ref_grid_tri(ref_grid);
  ncell += ref_cell_n(ref_cell);
  size += ref_cell_n(ref_cell) * (1 + ref_cell_node_per(ref_cell));

  ref_cell = ref_grid_qua(ref_grid);
  ncell += ref_cell_n(ref_cell);
  size += ref_cell_n(ref_cell) * (1 + ref_cell_node_per(ref_cell));

  each_ref_grid_ref_cell(ref_grid, group, ref_cell) {
    ncell += ref_cell_n(ref_cell);
    size += ref_cell_n(ref_cell) * (1 + ref_cell_node_per(ref_cell));
  }

  fprintf(file, "CELLS %d %d\n", ncell, size);

  ref_cell = ref_grid_tri(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    fprintf(file, " %d", node_per);
    for (node = 0; node < node_per; node++)
      fprintf(file, " %d", o2n[nodes[node]]);
    fprintf(file, "\n");
  }

  ref_cell = ref_grid_qua(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    fprintf(file, " %d", node_per);
    for (node = 0; node < node_per; node++)
      fprintf(file, " %d", o2n[nodes[node]]);
    fprintf(file, "\n");
  }

  each_ref_grid_ref_cell(ref_grid, group, ref_cell) {
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

REF_STATUS ref_export_tec(REF_GRID ref_grid, const char *filename) {
  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine geometry file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

  RSS(ref_export_tec_edge_zone(ref_grid, file), "edge");
  RSS(ref_export_tec_surf_zone(ref_grid, file), "surf");
  RSS(ref_export_tec_vol_zone(ref_grid, file), "vol");

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

REF_STATUS ref_export_tec_edge_zone(REF_GRID ref_grid, FILE *file) {
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

  ref_cell = ref_grid_edg(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes)
      RSS(ref_dict_store(ref_dict, nodes[2], REF_EMPTY), "mark tri");

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

    ref_cell = ref_grid_edg(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (boundary_tag == nodes[2]) {
        for (node = 0; node < 2; node++) {
          fprintf(file, " %d", g2l[nodes[node]] + 1);
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

REF_STATUS ref_export_tec_surf_zone(REF_GRID ref_grid, FILE *file) {
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
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes)
      RSS(ref_dict_store(ref_dict, nodes[3], REF_EMPTY), "mark tri");

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes)
      RSS(ref_dict_store(ref_dict, nodes[4], REF_EMPTY), "mark qua");

  each_ref_dict_key(ref_dict, boundary_index, boundary_tag) {
    RSS(ref_grid_tri_qua_id_nodes(ref_grid, boundary_tag, &nnode, &nface, &g2l,
                                  &l2g),
        "extract this boundary");

    fprintf(file,
            "zone t=\"surf%d\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
            boundary_tag, nnode, nface, "point", "fequadrilateral");

    for (node = 0; node < nnode; node++)
      fprintf(file, " %.16e %.16e %.16e\n",
              ref_node_xyz(ref_node, 0, l2g[node]),
              ref_node_xyz(ref_node, 1, l2g[node]),
              ref_node_xyz(ref_node, 2, l2g[node]));

    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell,
                                        nodes) if (boundary_tag == nodes[3]) {
      nodes[3] = nodes[2];
      for (node = 0; node < 4; node++) {
        fprintf(file, " %d", g2l[nodes[node]] + 1);
      }
      fprintf(file, "\n");
    }

    ref_cell = ref_grid_qua(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell,
                                        nodes) if (boundary_tag == nodes[4]) {
      for (node = 0; node < 4; node++)
        fprintf(file, " %d", g2l[nodes[node]] + 1);
      fprintf(file, "\n");
    }

    ref_free(l2g);
    ref_free(g2l);
  }

  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_vol_zone(REF_GRID ref_grid, FILE *file) {
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

  each_ref_grid_ref_cell(ref_grid, group,
                         ref_cell) if (ref_cell_n(ref_cell) > 0) {
    node_per = ref_cell_node_per(ref_cell);

    ref_malloc(o2n, ref_node_max(ref_node), REF_INT);

    for (node = 0; node < ref_node_max(ref_node); node++) o2n[node] = REF_EMPTY;

    /* mark nodes needed by this element type */
    nnode = 0;
    each_ref_cell_valid_cell_with_nodes(
        ref_cell, cell, nodes) for (node = 0; node < node_per;
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

    fprintf(
        file,
        "zone t=\"e%d\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
        node_per, nnode, ref_cell_n(ref_cell), "point", "febrick");

    for (node = 0; node < nnode; node++)
      fprintf(file, " %.16e %.16e %.16e\n",
              ref_node_xyz(ref_node, 0, n2o[node]),
              ref_node_xyz(ref_node, 1, n2o[node]),
              ref_node_xyz(ref_node, 2, n2o[node]));

    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      switch (ref_cell_node_per(ref_cell)) {
        case 4:
          TEC_BRICK_TET(brick, nodes);
          break;
        case 5:
          TEC_BRICK_PYR(brick, nodes);
          break;
        case 6:
          TEC_BRICK_PRI(brick, nodes);
          break;
        case 8:
          TEC_BRICK_HEX(brick, nodes);
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

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_int(REF_GRID ref_grid, REF_INT *scalar,
                              const char *filename) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT brick[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT ncell;
  REF_INT group;

  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine scalar file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\" \"s\"\n");

  ncell = 0;
  each_ref_grid_ref_cell(ref_grid, group, ref_cell) ncell +=
      ref_cell_n(ref_cell);

  fprintf(
      file,
      "zone t=\"scalar\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      ref_node_n(ref_node), ncell, "point", "febrick");

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  for (node = 0; node < ref_node_n(ref_node); node++)
    fprintf(file, " %.16e %.16e %.16e %d\n",
            ref_node_xyz(ref_node, 0, n2o[node]),
            ref_node_xyz(ref_node, 1, n2o[node]),
            ref_node_xyz(ref_node, 2, n2o[node]), scalar[n2o[node]]);

  each_ref_grid_ref_cell(ref_grid, group, ref_cell)
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    switch (ref_cell_node_per(ref_cell)) {
      case 4:
        TEC_BRICK_TET(brick, nodes);
        break;
      case 5:
        TEC_BRICK_PYR(brick, nodes);
        break;
      case 6:
        TEC_BRICK_PRI(brick, nodes);
        break;
      case 8:
        TEC_BRICK_HEX(brick, nodes);
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

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_dbl(REF_GRID ref_grid, REF_INT ldim, REF_DBL *scalar,
                              const char *filename) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT brick[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT ncell;
  REF_INT group;
  REF_INT i;

  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine scalar file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"");
  for (i = 0; i < ldim; i++) {
    fprintf(file, " \"s%d\"", i + 1);
  }
  fprintf(file, "\n");

  ncell = 0;
  each_ref_grid_ref_cell(ref_grid, group, ref_cell) ncell +=
      ref_cell_n(ref_cell);

  fprintf(
      file,
      "zone t=\"scalar\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      ref_node_n(ref_node), ncell, "point", "febrick");

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  for (node = 0; node < ref_node_n(ref_node); node++) {
    fprintf(file, " %.16e %.16e %.16e", ref_node_xyz(ref_node, 0, n2o[node]),
            ref_node_xyz(ref_node, 1, n2o[node]),
            ref_node_xyz(ref_node, 2, n2o[node]));
    for (i = 0; i < ldim; i++) {
      fprintf(file, " %.16e", scalar[i + ldim * node]);
    }
    fprintf(file, "\n");
  }

  each_ref_grid_ref_cell(ref_grid, group, ref_cell)
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    switch (ref_cell_node_per(ref_cell)) {
      case 4:
        TEC_BRICK_TET(brick, nodes);
        break;
      case 5:
        TEC_BRICK_PYR(brick, nodes);
        break;
      case 6:
        TEC_BRICK_PRI(brick, nodes);
        break;
      case 8:
        TEC_BRICK_HEX(brick, nodes);
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

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_tec_part(REF_GRID ref_grid, const char *root_filename) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  char viz_file[256];

  sprintf(viz_file, "%s_n%d_p%d.tec", root_filename,
          ref_mpi_n(ref_grid_mpi(ref_grid)),
          ref_mpi_rank(ref_grid_mpi(ref_grid)));

  RSS(ref_export_tec_int(ref_grid, ref_node->part, viz_file),
      "viz parts as scalar");

  return REF_SUCCESS;
}

REF_STATUS ref_export_metric_xyzdirlen(REF_GRID ref_grid,
                                       const char *filename) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_DBL m[6], d[12];
  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  for (node = 0; node < ref_node_n(ref_node); node++) {
    RSS(ref_node_metric_get(ref_node, n2o[node], m), "get");
    RSS(ref_matrix_diag_m(m, d), "diag");
    RSS(ref_matrix_ascending_eig(d), "sort eig");

    fprintf(file,
            " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e "
            "%.16e %.16e %.16e %.16e %.16e\n",
            ref_node_xyz(ref_node, 0, n2o[node]),
            ref_node_xyz(ref_node, 1, n2o[node]),
            ref_node_xyz(ref_node, 2, n2o[node]), d[3], d[4], d[5], d[6], d[7],
            d[8], d[9], d[10], d[11], 1.0 / sqrt(d[0]), 1.0 / sqrt(d[1]),
            1.0 / sqrt(d[2]));
  }

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}
REF_STATUS ref_export_tec_metric_axis(REF_GRID ref_grid,
                                      const char *root_filename) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT brick[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;
  REF_INT ncell;
  REF_INT group;
  REF_DBL m[6], d[12];
  REF_DBL dx, dy, dz;
  FILE *file;
  char viz_file[256];
  REF_INT e;

  for (e = 0; e < 3; e++) {
    sprintf(viz_file, "%s_n%d_p%d_ax%d.tec", root_filename,
            ref_mpi_n(ref_grid_mpi(ref_grid)),
            ref_mpi_rank(ref_grid_mpi(ref_grid)), e);

    file = fopen(viz_file, "w");
    if (NULL == (void *)file) printf("unable to open %s\n", viz_file);
    RNS(file, "unable to open file");

    fprintf(file, "title=\"tecplot refine metric axes\"\n");
    fprintf(file, "variables = \"x\" \"y\" \"z\" \"u\" \"v\" \"w\"\n");

    ncell = 0;
    each_ref_grid_ref_cell(ref_grid, group, ref_cell) ncell +=
        ref_cell_n(ref_cell);

    fprintf(file,
            "zone t=\"scalar\", nodes=%d, elements=%d, datapacking=%s, "
            "zonetype=%s\n",
            ref_node_n(ref_node), ncell, "point", "febrick");

    RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

    for (node = 0; node < ref_node_n(ref_node); node++) {
      RSS(ref_node_metric_get(ref_node, n2o[node], m), "get");
      RSS(ref_matrix_diag_m(m, d), "diag");
      RSS(ref_matrix_ascending_eig(d), "sort eig");
      dx = d[3 + 3 * e] / sqrt(d[e]);
      dy = d[4 + 3 * e] / sqrt(d[e]);
      dz = d[5 + 3 * e] / sqrt(d[e]);
      fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e\n",
              ref_node_xyz(ref_node, 0, n2o[node]),
              ref_node_xyz(ref_node, 1, n2o[node]),
              ref_node_xyz(ref_node, 2, n2o[node]), dx, dy, dz);
    }

    each_ref_grid_ref_cell(ref_grid, group, ref_cell)
        each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      switch (ref_cell_node_per(ref_cell)) {
        case 4:
          TEC_BRICK_TET(brick, nodes);
          break;
        case 5:
          TEC_BRICK_PYR(brick, nodes);
          break;
        case 6:
          TEC_BRICK_PRI(brick, nodes);
          break;
        case 8:
          TEC_BRICK_HEX(brick, nodes);
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

    fclose(file);
  } /* each eigenpair */

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
  REF_DBL x, y, z;
  REF_DBL ex, ey;
  FILE *file;
  char viz_file[256];
  REF_INT i, n = 36;
  REF_INT e0, e1, eb;
  REF_DBL best_y;
  REF_DBL dt = ref_math_in_radians(360.0 / (REF_DBL)n);
  REF_DBL scale = 0.5; /* so the ellipses touch for an ideal grid */

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

  fprintf(
      file,
      "zone t=\"scalar\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      1 * ncell, 1 * ncell, "point", "felineseg");

  for (node = 0; node < nnode; node++) {
    RSS(ref_node_metric_get(ref_node, n2o[node], m), "get");
    RSS(ref_matrix_diag_m(m, d), "diag");
    RSS(ref_matrix_ascending_eig(d), "sort eig");
    eb = REF_EMPTY;
    best_y = -1.0;
    for (e0 = 0; e0 < 3; e0++) {
      if (ABS(d[4 + 3 * e0]) > best_y) {
        eb = e0;
        best_y = ABS(d[4 + 3 * e0]);
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
      z = d[5 + 3 * e0] * ex + d[5 + 3 * e1] * ey;
      fprintf(file, " %.16e %.16e %.16e\n",
              ref_node_xyz(ref_node, 0, n2o[node]) + x,
              ref_node_xyz(ref_node, 1, n2o[node]) + y,
              ref_node_xyz(ref_node, 2, n2o[node]) + z);
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

  RSS(ref_export_tec_surf_zone(ref_grid, file), "ellipse surf");

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

  fprintf(
      file,
      "zone t=\"scalar\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
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

REF_STATUS ref_export_tec_ratio(REF_GRID ref_grid, const char *root_filename) {
  REF_EDGE ref_edge;
  char viz_file[256];

  sprintf(viz_file, "%s_n%d_p%d.tec", root_filename,
          ref_mpi_n(ref_grid_mpi(ref_grid)),
          ref_mpi_rank(ref_grid_mpi(ref_grid)));

  RSS(ref_edge_create(&ref_edge, ref_grid), "make edge");

  RSS(ref_edge_tec_ratio(ref_edge, viz_file), "viz parts as scalar");

  RSS(ref_edge_free(ref_edge), "free edge");

  return REF_SUCCESS;
}

REF_STATUS ref_export_poly(REF_GRID ref_grid, const char *filename) {
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
      RUS(REF_EMPTY, largest_triangle, "no largest triangle");
      RSS(ref_cell_nodes(ref_cell, largest_triangle, nodes), "tri nodes");
      RSS(ref_node_tri_normal(ref_node, nodes, normal), "normal");
      RSS(ref_math_normalize(normal), "norm");
      offset = 1.0e-4 * sqrt(area);
      if (offset > 1.0e-12) {
        nhole += 1;
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
      RUS(REF_EMPTY, largest_triangle, "no largest triangle");
      RSS(ref_cell_nodes(ref_cell, largest_triangle, nodes), "tri nodes");
      RSS(ref_node_tri_normal(ref_node, nodes, normal), "normal");
      RSS(ref_math_normalize(normal), "norm");
      offset = 1.0e-4 * sqrt(area);
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

  /* Part 4 - region attributes list */
  fprintf(file, "0\n"); /* no regions */

  ref_free(n2o);
  ref_free(o2n);

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_smesh(REF_GRID ref_grid, const char *filename) {
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

REF_STATUS ref_export_fgrid(REF_GRID ref_grid, const char *filename) {
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
REF_STATUS ref_export_su2(REF_GRID ref_grid, const char *filename) {
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

  fprintf(file, "NDIME= 3\n");

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  nnode = ref_node_n(ref_node);

  fprintf(file, "NPOIN= %d\n", nnode);

  for (node = 0; node < ref_node_n(ref_node); node++)
    fprintf(file, " %.16e %.16e %.16e\n", ref_node_xyz(ref_node, 0, n2o[node]),
            ref_node_xyz(ref_node, 1, n2o[node]),
            ref_node_xyz(ref_node, 2, n2o[node]));

  ntet = ref_cell_n(ref_grid_tet(ref_grid));
  npyr = ref_cell_n(ref_grid_pyr(ref_grid));
  npri = ref_cell_n(ref_grid_pri(ref_grid));
  nhex = ref_cell_n(ref_grid_hex(ref_grid));

  fprintf(file, "NELEM= %d\n", ntet + npyr + npri + nhex);

  /* uses VTK windings */
  each_ref_grid_ref_cell(ref_grid, group, ref_cell) {
    node_per = ref_cell_node_per(ref_cell);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      switch (ref_cell_node_per(ref_cell)) {
        case 4:
          fprintf(file, "10");
          break;
        case 5:
          fprintf(file, "14");
          VTK_PYRAMID_ORDER(nodes);
          break;
        case 6:
          fprintf(file, "13");
          VTK_WEDGE_ORDER(nodes);
          break;
        case 8:
          fprintf(file, "12");
          break;
        default:
          RSS(REF_IMPLEMENT, "wrong nodes per cell");
          break;
      }
      for (node = 0; node < node_per; node++)
        fprintf(file, " %d", o2n[nodes[node]]);
      fprintf(file, "\n");
    }
  }

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
  for (faceid = min_faceid; faceid <= max_faceid; faceid++)
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[node_per] == faceid) {
        fprintf(file, " %d", nodes[3]);
        fprintf(file, "\n");
      }
    }

  ref_cell = ref_grid_qua(ref_grid);
  node_per = ref_cell_node_per(ref_cell);
  for (faceid = min_faceid; faceid <= max_faceid; faceid++)
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (nodes[node_per] == faceid) {
        fprintf(file, " %d", nodes[4]);
        fprintf(file, "\n");
      }
    }

  each_ref_grid_ref_cell(ref_grid, group, ref_cell) {
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

  each_ref_grid_ref_cell(ref_grid, group, ref_cell) {
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

REF_STATUS ref_export_cogsg(REF_GRID ref_grid, const char *filename_cogsg) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node, cell, ixyz;
  REF_INT nnode, nbn, ntri;
  REF_INT *o2n, *n2o;
  REF_CELL ref_cell;
  REF_DICT ref_dict;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  FILE *file;
  REF_INT boundary_index, boundary_tag;

  int fortran_record_size;
  REF_INT i;
  REF_DBL t;
  size_t end_of_string;
  char filename_bc[1024];
  char *period_bc;
  end_of_string = strlen(filename_cogsg);

  RAS(end_of_string > 6, "filename too short");
  RAS(end_of_string < 1023, "filename too long");
  REIS(0, strcmp(&filename_cogsg[end_of_string - 6], ".cogsg"),
       "filename must end in .cogsg");

  sprintf(filename_bc, "%s", filename_cogsg);
  period_bc = strrchr(filename_bc, '.');
  sprintf(period_bc, ".bc");

  REIS(0, ref_cell_n(ref_grid_qua(ref_grid)), "no quad support");
  REIS(0, ref_cell_n(ref_grid_pyr(ref_grid)), "no pyramid support");
  REIS(0, ref_cell_n(ref_grid_pri(ref_grid)), "no prism support");
  REIS(0, ref_cell_n(ref_grid_hex(ref_grid)), "no hex support");

  ref_malloc_init(o2n, ref_node_max(ref_node), REF_INT, REF_EMPTY);
  ref_malloc(n2o, ref_node_n(ref_node), REF_INT);

  nnode = 0;

  ref_cell = ref_grid_tri(ref_grid);

  each_ref_cell_valid_cell_with_nodes(
      ref_cell, cell, nodes) for (node = 0; node < ref_cell_node_per(ref_cell);
                                  node++) if (REF_EMPTY == o2n[nodes[node]]) {
    o2n[nodes[node]] = nnode;
    nnode++;
  }
  nbn = nnode;

  each_ref_node_valid_node(ref_node, node) if (REF_EMPTY == o2n[node]) {
    o2n[node] = nnode;
    nnode++;
  }

  REIS(nnode, ref_node_n(ref_node), "nnode miscount");

  each_ref_node_valid_node(ref_node, node) n2o[o2n[node]] = node;

  RSS(ref_dict_create(&ref_dict), "create dict");

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes)
      RSS(ref_dict_store(ref_dict, nodes[3], REF_EMPTY), "mark tri");

  /*
  2)  project.bc: patch/surface-triangle file (ASCII)
       integer fapch(mbf),fnode(mbf,3)
       character*1 text(80)
       open(unit=12,file='project.bc',form='formatted')

       read(12,*)nbf,nbc,npatch,igrid
       read(12,900)text
       do if=1,nbf
        read(12,*)jf,fapch(if),fnode(if,in),in=1,3)
       enddo
  900  format(80a1)

        nbf = number of boundary triangular faces
        nbc = number of surface grid nodes along patch boundaries (curves)
     npatch = number of surface patches
      igrid = 1 for inviscid grids; 2 for viscous grids
       text = text line
         jf = triangle index
  fapch(if) = surface patch index containing surface triangle "if"
fnode(if,in) = node "in" of triangle "if"

Note: triangle connectivities are according to the right-hand rule with
      the outward normals pointing into the computational domain.
   */

  file = fopen(filename_bc, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename_bc);
  RNS(file, "unable to open file");

  fprintf(file, "%d %d %d %d\n", ref_cell_n(ref_cell), 0, ref_dict_n(ref_dict),
          1);
  fprintf(file,
          "exported by ref_export_cogsg x         x         x         x        "
          " x        80\n");

  ntri = 0;

  each_ref_dict_key(ref_dict, boundary_index, boundary_tag)
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell,
                                          nodes) if (boundary_tag == nodes[3]) {
    fprintf(file, " %d %d", ntri + 1, boundary_tag);
    ntri++;
    for (node = 0; node < 3; node++) {
      fprintf(file, " %d", o2n[nodes[node]] + 1);
    }
    fprintf(file, "\n");
  }

  fclose(file);

  RSS(ref_dict_free(ref_dict), "free");

  /*
  3)  project.cogsg: x,y,z coordinates of the grid nodes and
                    tetrahedral node connectivity
      real(8)  crd(mp,3),t
      integer int(mc,4)
       open(9,file='project.cogrd',form='unformatted',iostat=ios,
     &      err=555,status='old')

      read(9)inew,nc,npo,nbn,npv,nev,t,
     &       ((int(ie,in),ie=1,nc),in=1,4)
      read(9)((crd(ip,id),ip=1,npo),id=1,3)
      read(9)0
      read(9)0

  where

         inew = a dummy variable (integer) should be -1
           nc = number of tetrahedral cells
          npo = total number of grid nodes (including nbn)
          nbn = number of grid nodes on the boundaries (including nbc)
          npv = number of grid points in the viscous layers
                (=0 for Euler grids)
          ncv = number of cells in the viscous layers
                (=0 for Euler grids)
            t = a dummy variable (real - double)
   int(ie,in) = tetradhedral cell connectivity
                (node "in" of cell "ie")
   crd(ip,id) = x, y, and z coordinates of node "ip"

  Note 1: the first "nbn" coordinates listed in this file are those of
         the boundary nodes.

  Note 2:  tetrahedral cell connectivities are given according to the
           right-hand rule (3 nodes of the base in the counter-clockwise
           direction followed by the 4th node.)
   */

  file = fopen(filename_cogsg, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename_cogsg);
  RNS(file, "unable to open file");

  ref_cell = ref_grid_tet(ref_grid);

  fortran_record_size = 4 * 6 + 8 * 1 + 4 * 4 * ref_cell_n(ref_cell);

  i = fortran_record_size;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "record");

  i = -1;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "inew");
  i = ref_cell_n(ref_cell);
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "nc");
  i = nnode;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "npo");
  i = nbn;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "nbn");
  i = 0;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "npv");
  i = 0;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "nev");
  t = 0.0;
  SWAP_DBL(t);
  REIS(1, fwrite(&t, sizeof(REF_DBL), 1, file), "t");

  for (node = 0; node < ref_cell_node_per(ref_cell); node++)
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      i = o2n[nodes[node]] + 1;
      SWAP_INT(i);
      REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "int");
    }

  i = fortran_record_size;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "record");

  fortran_record_size = 8 * 3 * nnode;

  i = fortran_record_size;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "record");
  for (ixyz = 0; ixyz < 3; ixyz++)
    for (node = 0; node < nnode; node++) {
      t = ref_node_xyz(ref_node, ixyz, n2o[node]);
      SWAP_DBL(t);
      REIS(1, fwrite(&t, sizeof(REF_DBL), 1, file), "crd");
    }
  i = fortran_record_size;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "record");

  i = 4;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "record");
  i = 0;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "nelec");
  i = 4;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "record");

  i = 4;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "record");
  i = 0;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "npoic");
  i = 4;
  SWAP_INT(i);
  REIS(1, fwrite(&i, sizeof(REF_INT), 1, file), "record");

  fclose(file);

  ref_free(o2n);
  ref_free(n2o);

  return REF_SUCCESS;
}

REF_STATUS ref_export_c(REF_GRID ref_grid, const char *filename) {
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

REF_STATUS ref_export_pdf(REF_GRID ref_grid, const char *pdf_filename) {
  char temp_filename[] = "ref_export_temp_for_pdf.eps";
  char command[1024];
  RSS(ref_export_eps(ref_grid, temp_filename), "temp eps");
  sprintf(command, "epstopdf %s -o=%s", temp_filename, pdf_filename);
  REIS(0, system(command), "epstopdf failed");
  REIS(0, remove(temp_filename), "temp clean up");

  return REF_SUCCESS;
}

REF_STATUS ref_export_eps(REF_GRID ref_grid, const char *filename) {
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

REF_STATUS ref_export_html(REF_GRID ref_grid, const char *filename) {
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
      printf("next_position outside int32 limits %d %d\n", -2147483647,
             2147483647);
      RSS(REF_INVALID, "meshb version does not support file size");
    }
    one_word = (int32_t)next_position;
    REIS(1, fwrite(&one_word, sizeof(one_word), 1, file), "write next pos");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_export_meshb(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_CELL ref_cell;
  REF_INT *o2n, *n2o;
  REF_INT code, version, dim;
  REF_FILEPOS next_position;
  REF_INT keyword_code, header_size;
  REF_INT node;
  REF_INT min_faceid, max_faceid, node_per, faceid, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT type, id, i, size_bytes;
  REF_INT geom;
  int ngeom;

  if (10000000 < ref_node_n(ref_node)) {
    version = 3;
    header_size = 4 + 8 + 4;
  } else {
    version = 2;
    header_size = 4 + 4 + 4;
  }

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  code = 1;
  REIS(1, fwrite(&code, sizeof(int), 1, file), "code");
  REIS(1, fwrite(&version, sizeof(int), 1, file), "version");
  next_position = header_size + ftell(file);
  keyword_code = 3;
  REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "dim code");
  RSS(ref_export_meshb_next_position(file, version, next_position), "next pos");
  dim = 3;
  if (ref_grid_twod(ref_grid)) dim = 2;
  REIS(1, fwrite(&dim, sizeof(int), 1, file), "dim");
  REIS(next_position, ftell(file), "dim inconsistent");

  if (ref_node_n(ref_node) > 0) {
    next_position =
        (REF_FILEPOS)header_size +
        (REF_FILEPOS)ref_node_n(ref_node) * (REF_FILEPOS)(dim * 8 + 4) +
        ftell(file);
    keyword_code = 4;
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "vertex version code");
    RSS(ref_export_meshb_next_position(file, version, next_position), "next p");
    REIS(1, fwrite(&(ref_node_n(ref_node)), sizeof(int), 1, file), "nnode");
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
      id = REF_EXPORT_MESHB_VERTEX_ID;
      REIS(1, fwrite(&(id), sizeof(int), 1, file), "id");
    }
    REIS(next_position, ftell(file), "vertex inconsistent");
  }

  ref_cell = ref_grid_edg(ref_grid);
  keyword_code = 5;
  if (ref_cell_n(ref_cell) > 0) {
    node_per = ref_cell_node_per(ref_cell);
    next_position =
        (REF_FILEPOS)header_size + ftell(file) +
        (REF_FILEPOS)ref_cell_n(ref_cell) * (REF_FILEPOS)(4 * (node_per + 1));
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "vertex version code");
    RSS(ref_export_meshb_next_position(file, version, next_position), "next p");
    REIS(1, fwrite(&(ref_cell_n(ref_cell)), sizeof(int), 1, file), "nnode");
    RSS(ref_export_edgeid_range(ref_grid, &min_faceid, &max_faceid), "range");
    for (faceid = min_faceid; faceid <= max_faceid; faceid++)
      each_ref_cell_valid_cell_with_nodes(
          ref_cell, cell, nodes) if (nodes[node_per] == faceid) {
        for (node = 0; node < node_per; node++) {
          nodes[node] = o2n[nodes[node]] + 1;
          REIS(1, fwrite(&(nodes[node]), sizeof(REF_INT), 1, file), "ele");
        }
        REIS(1, fwrite(&(nodes[2]), sizeof(REF_INT), 1, file), "ele id");
      }
    REIS(next_position, ftell(file), "edge inconsistent");
  }

  ref_cell = ref_grid_tri(ref_grid);
  keyword_code = 6;
  if (ref_cell_n(ref_cell) > 0) {
    node_per = ref_cell_node_per(ref_cell);
    next_position =
        (REF_FILEPOS)header_size + ftell(file) +
        (REF_FILEPOS)ref_cell_n(ref_cell) * (REF_FILEPOS)(4 * (node_per + 1));
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "vertex version code");
    RSS(ref_export_meshb_next_position(file, version, next_position), "next p");
    REIS(1, fwrite(&(ref_cell_n(ref_cell)), sizeof(int), 1, file), "nnode");
    RSS(ref_export_faceid_range(ref_grid, &min_faceid, &max_faceid), "range");

    for (faceid = min_faceid; faceid <= max_faceid; faceid++)
      each_ref_cell_valid_cell_with_nodes(
          ref_cell, cell, nodes) if (nodes[node_per] == faceid) {
        for (node = 0; node < node_per; node++) {
          nodes[node] = o2n[nodes[node]] + 1;
          REIS(1, fwrite(&(nodes[node]), sizeof(REF_INT), 1, file), "ele");
        }
        REIS(1, fwrite(&(nodes[3]), sizeof(REF_INT), 1, file), "ele id");
      }
    REIS(next_position, ftell(file), "tri inconsistent");
  }

  ref_cell = ref_grid_qua(ref_grid);
  keyword_code = 7;
  if (ref_cell_n(ref_cell) > 0) {
    node_per = ref_cell_node_per(ref_cell);
    next_position =
        (REF_FILEPOS)header_size + ftell(file) +
        (REF_FILEPOS)ref_cell_n(ref_cell) * (REF_FILEPOS)(5 * (node_per + 1));
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "vertex version code");
    RSS(ref_export_meshb_next_position(file, version, next_position), "next p");
    REIS(1, fwrite(&(ref_cell_n(ref_cell)), sizeof(int), 1, file), "nnode");
    RSS(ref_export_faceid_range(ref_grid, &min_faceid, &max_faceid), "range");

    for (faceid = min_faceid; faceid <= max_faceid; faceid++)
      each_ref_cell_valid_cell_with_nodes(
          ref_cell, cell, nodes) if (nodes[node_per] == faceid) {
        for (node = 0; node < node_per; node++) {
          nodes[node] = o2n[nodes[node]] + 1;
          REIS(1, fwrite(&(nodes[node]), sizeof(REF_INT), 1, file), "ele");
        }
        REIS(1, fwrite(&(nodes[3]), sizeof(REF_INT), 1, file), "ele id");
      }
    REIS(next_position, ftell(file), "qua inconsistent");
  }

  ref_cell = ref_grid_tet(ref_grid);
  if (ref_cell_n(ref_cell) > 0) {
    node_per = ref_cell_node_per(ref_cell);
    next_position =
        (REF_FILEPOS)header_size +
        (REF_FILEPOS)ref_cell_n(ref_cell) * (REF_FILEPOS)(4 * (node_per + 1)) +
        ftell(file);
    keyword_code = 8;
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "vertex version code");
    RSS(ref_export_meshb_next_position(file, version, next_position), "next p");
    REIS(1, fwrite(&(ref_cell_n(ref_cell)), sizeof(int), 1, file), "nnode");
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      for (node = 0; node < node_per; node++) {
        nodes[node] = o2n[nodes[node]] + 1;
        REIS(1, fwrite(&(nodes[node]), sizeof(int), 1, file), "cell");
      }
      id = REF_EXPORT_MESHB_3D_ID;
      REIS(1, fwrite(&(id), sizeof(int), 1, file), "tet id");
    }
    REIS(next_position, ftell(file), "tet inconsistent");
  }

  each_ref_type(ref_geom, type) {
    ngeom = 0;
    each_ref_geom_of(ref_geom, type, geom) ngeom++;
    if (ngeom > 0) {
      keyword_code = 40 + type; /* GmfVerticesOnGeometricVertices */
      next_position = (REF_FILEPOS)header_size +
                      (REF_FILEPOS)ngeom *
                          (REF_FILEPOS)(4 * 2 + 8 * type + (0 < type ? 8 : 0)) +
                      ftell(file);
      REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "keyword");
      RSS(ref_export_meshb_next_position(file, version, next_position), "next");
      REIS(1, fwrite(&(ngeom), sizeof(int), 1, file), "n");
      each_ref_geom_of(ref_geom, type, geom) {
        double filler = 0.0;
        node = o2n[ref_geom_node(ref_geom, geom)] + 1;
        id = ref_geom_id(ref_geom, geom);
        REIS(1, fwrite(&(node), sizeof(int), 1, file), "node");
        REIS(1, fwrite(&(id), sizeof(int), 1, file), "id");
        for (i = 0; i < type; i++)
          REIS(1,
               fwrite(&(ref_geom_param(ref_geom, i, geom)), sizeof(double), 1,
                      file),
               "id");
        if (0 < type) REIS(1, fwrite(&(filler), sizeof(double), 1, file), "id");
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
    REIS(1, fwrite(&size_bytes, sizeof(int), 1, file), "n");
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

REF_STATUS ref_export_twod_msh(REF_GRID ref_grid, const char *filename) {
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
  fprintf(f, "\nTriangles\n%d\n", ref_cell_n(ref_cell) / 2);
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

REF_STATUS ref_export_metric2d(REF_GRID ref_grid, const char *filename) {
  FILE *f;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nnode;
  REF_DBL m[6];

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

  fprintf(f, "%d %d\n", nnode, 3);
  for (node = 0; node < nnode; node++) {
    RSS(ref_node_metric_get(ref_node, n2o[node], m), "get");
    fprintf(f, "%.16E %.16E %.16E \n", m[0], m[2], m[5]);
  }

  ref_free(n2o);
  ref_free(o2n);

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_export_twod_sol(REF_GRID ref_grid, const char *filename) {
  FILE *f;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_INT *o2n, *n2o;
  REF_INT nnode;
  REF_DBL m[6];

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

  fprintf(f, "MeshVersionFormatted 2\n\n");
  fprintf(f, "Dimension 2\n\n");
  fprintf(f, "SolAtVertices\n%d\n1 3\n", nnode);

  for (node = 0; node < nnode; node++) {
    RSS(ref_node_metric_get(ref_node, n2o[node], m), "get");
    fprintf(f, "%.16E %.16E %.16E \n", m[0], m[2], m[5]);
  }

  ref_free(n2o);
  ref_free(o2n);

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_export_plt(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  int one = 1;
  int filetype = 0;
  int ascii[8];
  int numvar = 3;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  REIS(8, fwrite(&"#!TDV112", sizeof(char), 8, file), "header");
  REIS(1, fwrite(&one, sizeof(int), 1, file), "magic");
  REIS(1, fwrite(&filetype, sizeof(int), 1, file), "filetype");

  ascii[0] = (int)'f';
  ascii[1] = (int)'t';
  ascii[2] = 0;
  REIS(3, fwrite(&ascii, sizeof(int), 3, file), "title");

  REIS(1, fwrite(&numvar, sizeof(int), 1, file), "numvar");
  ascii[0] = (int)'x';
  ascii[1] = 0;
  REIS(2, fwrite(&ascii, sizeof(int), 2, file), "var");
  ascii[0] = (int)'y';
  ascii[1] = 0;
  REIS(2, fwrite(&ascii, sizeof(int), 2, file), "var");
  ascii[0] = (int)'z';
  ascii[1] = 0;
  REIS(2, fwrite(&ascii, sizeof(int), 2, file), "var");

  RSS(ref_export_plt_tet_zone(ref_grid, file), "plt tet zone");

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_export_plt_tet_zone(REF_GRID ref_grid, FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  int ascii[8];
  float zonemarker = 299.0;
  int parentzone = -1;
  int strandid = -1;
  double solutiontime = 0.0;
  int notused = -1;
  int zonetype = 4;    /*4=FETETRAHEDRON*/
  int datapacking = 0; /*0=Block, point does not work.*/
  int varloc = 0;      /*0 = Don't specify, all data is located at nodes*/
  int faceneighbors = 0;
  int numpts = ref_node_n(ref_node);
  int numelements = ref_cell_n(ref_cell);
  int celldim = 0;
  int aux = 0;
  float eohmarker = 357.0;
  int dataformat = 1;
  int passive = 0;
  int varsharing = 0;
  int connsharing = -1;
  float data;
  double mindata, maxdata;
  REF_INT *o2n, *n2o;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell, node, node_per, ixyz;
  int index;

  REIS(1, fwrite(&zonemarker, sizeof(float), 1, file), "zonemarker");

  ascii[0] = (int)'e';
  ascii[1] = (int)'4';
  ascii[2] = 0;
  REIS(3, fwrite(&ascii, sizeof(int), 3, file), "title");

  REIS(1, fwrite(&parentzone, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&strandid, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&solutiontime, sizeof(double), 1, file), "double");
  REIS(1, fwrite(&notused, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&zonetype, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&datapacking, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&varloc, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&faceneighbors, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&numpts, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&numelements, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&celldim, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&celldim, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&celldim, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&aux, sizeof(int), 1, file), "int");

  REIS(1, fwrite(&eohmarker, sizeof(float), 1, file), "eohmarker");
  REIS(1, fwrite(&zonemarker, sizeof(float), 1, file), "zonemarker");

  REIS(1, fwrite(&dataformat, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&dataformat, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&dataformat, sizeof(int), 1, file), "int");

  REIS(1, fwrite(&passive, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&varsharing, sizeof(int), 1, file), "int");
  REIS(1, fwrite(&connsharing, sizeof(int), 1, file), "int");

  RSS(ref_node_compact(ref_node, &o2n, &n2o), "compact");

  for (ixyz = 0; ixyz < 3; ixyz++) {
    mindata = ref_node_xyz(ref_node, ixyz, n2o[0]);
    maxdata = ref_node_xyz(ref_node, ixyz, n2o[0]);
    for (node = 1; node < ref_node_n(ref_node); node++) {
      mindata = MIN(mindata, ref_node_xyz(ref_node, ixyz, n2o[node]));
      maxdata = MAX(maxdata, ref_node_xyz(ref_node, ixyz, n2o[node]));
    }
    REIS(1, fwrite(&mindata, sizeof(double), 1, file), "mindata");
    REIS(1, fwrite(&maxdata, sizeof(double), 1, file), "maxdata");
  }

  for (ixyz = 0; ixyz < 3; ixyz++)
    for (node = 0; node < ref_node_n(ref_node); node++) {
      data = (float)ref_node_xyz(ref_node, ixyz, n2o[node]);
      REIS(1, fwrite(&data, sizeof(float), 1, file), "data");
    }

  node_per = ref_cell_node_per(ref_cell);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < node_per; node++) {
      index = o2n[nodes[node]];
      REIS(1, fwrite(&index, sizeof(int), 1, file), "index");
    }
  }

  ref_free(n2o);
  ref_free(o2n);

  return REF_SUCCESS;
}

REF_STATUS ref_export_plt_surf_zone(REF_GRID ref_grid, FILE *file) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  int ascii[8];
  float zonemarker = 299.0;
  int parentzone = -1;
  int strandid = -1;
  double solutiontime = 0.0;
  int notused = -1;
  int zonetype = 3;    /*3=FEQUADRILATERAL*/
  int datapacking = 0; /*0=Block, point does not work.*/
  int varloc = 0;      /*0 = Don't specify, all data is located at nodes*/
  int faceneighbors = 0;
  int numpts = ref_node_n(ref_node);
  int numelements = ref_cell_n(ref_cell);
  int celldim = 0;
  int aux = 0;
  float eohmarker = 357.0;
  int dataformat = 1;
  int passive = 0;
  int varsharing = 0;
  int connsharing = -1;
  float data;
  double mindata, maxdata;
  REF_INT *g2l, *l2g;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell, node, ixyz;
  int index;

  REF_DICT ref_dict;
  REF_INT boundary_tag, boundary_index;
  REF_INT nnode, nface;

  RSS(ref_dict_create(&ref_dict), "create dict");

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes)
      RSS(ref_dict_store(ref_dict, nodes[3], REF_EMPTY), "mark tri");

  each_ref_dict_key(ref_dict, boundary_index, boundary_tag) {
    RSS(ref_grid_tri_qua_id_nodes(ref_grid, boundary_tag, &nnode, &nface, &g2l,
                                  &l2g),
        "extract this boundary");

    REIS(1, fwrite(&zonemarker, sizeof(float), 1, file), "zonemarker");

    ascii[0] = (int)'s';
    ascii[1] = 0;
    REIS(2, fwrite(&ascii, sizeof(int), 2, file), "title");

    numpts = nnode;
    numelements = nface;

    REIS(1, fwrite(&parentzone, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&strandid, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&solutiontime, sizeof(double), 1, file), "double");
    REIS(1, fwrite(&notused, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&zonetype, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&datapacking, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&varloc, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&faceneighbors, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&numpts, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&numelements, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&celldim, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&celldim, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&celldim, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&aux, sizeof(int), 1, file), "int");

    REIS(1, fwrite(&eohmarker, sizeof(float), 1, file), "eohmarker");
    REIS(1, fwrite(&zonemarker, sizeof(float), 1, file), "zonemarker");

    REIS(1, fwrite(&dataformat, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&dataformat, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&dataformat, sizeof(int), 1, file), "int");

    REIS(1, fwrite(&passive, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&varsharing, sizeof(int), 1, file), "int");
    REIS(1, fwrite(&connsharing, sizeof(int), 1, file), "int");

    for (ixyz = 0; ixyz < 3; ixyz++) {
      mindata = ref_node_xyz(ref_node, ixyz, l2g[0]);
      maxdata = ref_node_xyz(ref_node, ixyz, l2g[0]);
      for (node = 1; node < nnode; node++) {
        mindata = MIN(mindata, ref_node_xyz(ref_node, ixyz, l2g[node]));
        maxdata = MAX(maxdata, ref_node_xyz(ref_node, ixyz, l2g[node]));
      }
      REIS(1, fwrite(&mindata, sizeof(double), 1, file), "mindata");
      REIS(1, fwrite(&maxdata, sizeof(double), 1, file), "maxdata");
    }

    for (ixyz = 0; ixyz < 3; ixyz++)
      for (node = 0; node < nnode; node++) {
        data = (float)ref_node_xyz(ref_node, ixyz, l2g[node]);
        REIS(1, fwrite(&data, sizeof(float), 1, file), "data");
      }

    ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell,
                                        nodes) if (boundary_tag == nodes[3]) {
      nodes[3] = nodes[2];
      for (node = 0; node < 4; node++) {
        index = g2l[nodes[node]];
        REIS(1, fwrite(&index, sizeof(int), 1, file), "index");
      }
    }

    ref_cell = ref_grid_qua(ref_grid);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell,
                                        nodes) if (boundary_tag == nodes[4]) {
      for (node = 0; node < 4; node++) {
        index = g2l[nodes[node]];
        REIS(1, fwrite(&index, sizeof(int), 1, file), "index");
      }
    }

    ref_free(l2g);
    ref_free(g2l);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_export_faceid_range(REF_GRID ref_grid, REF_INT *min_faceid,
                                   REF_INT *max_faceid) {
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  *min_faceid = REF_INT_MAX;
  *max_faceid = REF_INT_MIN;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    *min_faceid = MIN(*min_faceid, nodes[ref_cell_node_per(ref_cell)]);
    *max_faceid = MAX(*max_faceid, nodes[ref_cell_node_per(ref_cell)]);
  }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    *min_faceid = MIN(*min_faceid, nodes[ref_cell_node_per(ref_cell)]);
    *max_faceid = MAX(*max_faceid, nodes[ref_cell_node_per(ref_cell)]);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_export_edgeid_range(REF_GRID ref_grid, REF_INT *min_edgeid,
                                   REF_INT *max_edgeid) {
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  *min_edgeid = REF_INT_MAX;
  *max_edgeid = REF_INT_MIN;

  ref_cell = ref_grid_edg(ref_grid);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    *min_edgeid = MIN(*min_edgeid, nodes[ref_cell_node_per(ref_cell)]);
    *max_edgeid = MAX(*max_edgeid, nodes[ref_cell_node_per(ref_cell)]);
  }

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
    RSS(ref_export_plt(ref_grid, filename), "plt export failed");
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
  } else if (strcmp(&filename[end_of_string - 5], ".poly") == 0) {
    RSS(ref_export_poly(ref_grid, filename), "poly export failed");
  } else if (strcmp(&filename[end_of_string - 6], ".smesh") == 0) {
    RSS(ref_export_smesh(ref_grid, filename), "smesh export failed");
  } else if (strcmp(&filename[end_of_string - 6], ".fgrid") == 0) {
    RSS(ref_export_fgrid(ref_grid, filename), "fgrid export failed");
  } else if (strcmp(&filename[end_of_string - 6], ".cogsg") == 0) {
    RSS(ref_export_cogsg(ref_grid, filename), "cogsg export failed");
  } else if (strcmp(&filename[end_of_string - 5], ".html") == 0) {
    RSS(ref_export_html(ref_grid, filename), "html export failed");
  } else if (strcmp(&filename[end_of_string - 6], ".meshb") == 0) {
    RSS(ref_export_meshb(ref_grid, filename), "meshb export failed");
  } else if (strcmp(&filename[end_of_string - 4], ".msh") == 0) {
    RSS(ref_export_twod_msh(ref_grid, filename), "msh export failed");
  } else {
    printf("%s: %d: %s %s\n", __FILE__, __LINE__,
           "export file name extension unknown", filename);
    RSS(REF_FAILURE, "unknown file extension");
  }

  return REF_SUCCESS;
}
