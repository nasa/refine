
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

#include "ref_import.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_endian.h"
#include "ref_malloc.h"

#define VTK_LINE (3)
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

#define VTK_PYRAMID_TO_UGRID(vtk_nodes) \
  {                                     \
    REF_INT ugrid_nodes[5];             \
    ugrid_nodes[0] = (vtk_nodes)[1];    \
    ugrid_nodes[1] = (vtk_nodes)[0];    \
    ugrid_nodes[2] = (vtk_nodes)[4];    \
    ugrid_nodes[3] = (vtk_nodes)[2];    \
    ugrid_nodes[4] = (vtk_nodes)[3];    \
    (vtk_nodes)[0] = ugrid_nodes[0];    \
    (vtk_nodes)[1] = ugrid_nodes[1];    \
    (vtk_nodes)[2] = ugrid_nodes[2];    \
    (vtk_nodes)[3] = ugrid_nodes[3];    \
    (vtk_nodes)[4] = ugrid_nodes[4];    \
  }

/*
 /3-/0
5-+-2| UGRID
 \4-\1
 /4-/1
5-+-2| VTK
 \3-\0
*/
#define VTK_WEDGE_TO_UGRID(vtk_nodes) \
  {                                   \
    REF_INT ugrid_nodes[6];           \
    ugrid_nodes[0] = (vtk_nodes)[1];  \
    ugrid_nodes[1] = (vtk_nodes)[0];  \
    ugrid_nodes[2] = (vtk_nodes)[2];  \
    ugrid_nodes[3] = (vtk_nodes)[4];  \
    ugrid_nodes[4] = (vtk_nodes)[3];  \
    ugrid_nodes[5] = (vtk_nodes)[5];  \
    (vtk_nodes)[0] = ugrid_nodes[0];  \
    (vtk_nodes)[1] = ugrid_nodes[1];  \
    (vtk_nodes)[2] = ugrid_nodes[2];  \
    (vtk_nodes)[3] = ugrid_nodes[3];  \
    (vtk_nodes)[4] = ugrid_nodes[4];  \
    (vtk_nodes)[5] = ugrid_nodes[5];  \
  }

static REF_STATUS ref_import_fgrid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                   const char *filename) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  FILE *file;
  REF_INT nnode, ntri, ntet;
  REF_INT ixyz, node, new_node;
  REF_DBL xyz;
  REF_INT tri, new_tri;
  REF_INT nodes[4];
  REF_INT face_id;
  REF_INT cell, new_cell;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RES(1, fscanf(file, "%d", &nnode), "nnode");
  RES(1, fscanf(file, "%d", &ntri), "ntri");
  RES(1, fscanf(file, "%d", &ntet), "ntet");

  for (node = 0; node < nnode; node++) {
    RSS(ref_node_add(ref_node, node, &new_node), "new_node");
    RES(node, new_node, "node index");
  }

  RSS(ref_node_initialize_n_global(ref_node, nnode), "init glob");

  for (ixyz = 0; ixyz < 3; ixyz++)
    for (node = 0; node < nnode; node++) {
      RES(1, fscanf(file, "%lf", &xyz), "xyz");
      ref_node_xyz(ref_node, ixyz, node) = xyz;
    }

  ref_cell = ref_grid_tri(ref_grid);
  nodes[3] = REF_EMPTY;
  for (tri = 0; tri < ntri; tri++) {
    for (node = 0; node < 3; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tri");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_tri), "new tri");
    RES(tri, new_tri, "tri index");
  }

  ref_cell = ref_grid_tri(ref_grid);
  for (tri = 0; tri < ntri; tri++) {
    RES(1, fscanf(file, "%d", &face_id), "tri id");
    ref_cell_c2n(ref_cell, 3, tri) = face_id;
  }

  ref_cell = ref_grid_tet(ref_grid);
  for (cell = 0; cell < ntet; cell++) {
    for (node = 0; node < 4; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tet");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    nodes[3]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "new tet");
    RES(cell, new_cell, "tet index");
  }

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_import_ugrid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                   const char *filename) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_INT node, new_node;
  REF_DBL xyz[3];
  REF_INT tri, new_tri;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT qua, new_qua;
  REF_INT face_id;
  REF_INT cell, new_cell;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RES(1, fscanf(file, "%d", &nnode), "nnode");
  RES(1, fscanf(file, "%d", &ntri), "ntri");
  RES(1, fscanf(file, "%d", &nqua), "nqua");
  RES(1, fscanf(file, "%d", &ntet), "ntet");
  RES(1, fscanf(file, "%d", &npyr), "npyr");
  RES(1, fscanf(file, "%d", &npri), "npri");
  RES(1, fscanf(file, "%d", &nhex), "nhex");

  for (node = 0; node < nnode; node++) {
    RSS(ref_node_add(ref_node, node, &new_node), "new_node");
    RES(node, new_node, "node index");
    RES(1, fscanf(file, "%lf", &(xyz[0])), "x");
    RES(1, fscanf(file, "%lf", &(xyz[1])), "y");
    RES(1, fscanf(file, "%lf", &(xyz[2])), "z");
    ref_node_xyz(ref_node, 0, new_node) = xyz[0];
    ref_node_xyz(ref_node, 1, new_node) = xyz[1];
    ref_node_xyz(ref_node, 2, new_node) = xyz[2];
  }

  RSS(ref_node_initialize_n_global(ref_node, nnode), "init glob");

  ref_cell = ref_grid_tri(ref_grid);
  nodes[3] = REF_EMPTY;
  for (tri = 0; tri < ntri; tri++) {
    for (node = 0; node < 3; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tri");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_tri), "new tri");
    RES(tri, new_tri, "tri index");
  }

  ref_cell = ref_grid_qua(ref_grid);
  nodes[4] = REF_EMPTY;
  for (qua = 0; qua < nqua; qua++) {
    for (node = 0; node < 4; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "qua");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    nodes[3]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_qua), "new qua");
    RES(qua, new_qua, "qua index");
  }

  ref_cell = ref_grid_tri(ref_grid);
  for (tri = 0; tri < ntri; tri++) {
    RES(1, fscanf(file, "%d", &face_id), "tri id");
    ref_cell_c2n(ref_cell, 3, tri) = face_id;
  }

  ref_cell = ref_grid_qua(ref_grid);
  for (qua = 0; qua < nqua; qua++) {
    RES(1, fscanf(file, "%d", &face_id), "qua id");
    ref_cell_c2n(ref_cell, 4, qua) = face_id;
  }

  ref_cell = ref_grid_tet(ref_grid);
  for (cell = 0; cell < ntet; cell++) {
    for (node = 0; node < 4; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tet");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    nodes[3]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "new tet");
    RES(cell, new_cell, "tet index");
  }

  ref_cell = ref_grid_pyr(ref_grid);
  for (cell = 0; cell < npyr; cell++) {
    for (node = 0; node < 5; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "pyr");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    nodes[3]--;
    nodes[4]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "new pyr");
    RES(cell, new_cell, "pyr index");
  }

  ref_cell = ref_grid_pri(ref_grid);
  for (cell = 0; cell < npri; cell++) {
    for (node = 0; node < 6; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "pri");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    nodes[3]--;
    nodes[4]--;
    nodes[5]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "new pri");

    if (cell != new_cell) {
      printf("cell %d %d\n", cell, new_cell);
    }

    RES(cell, new_cell, "pri index");
  }

  ref_cell = ref_grid_hex(ref_grid);
  for (cell = 0; cell < nhex; cell++) {
    for (node = 0; node < 8; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "hex");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    nodes[3]--;
    nodes[4]--;
    nodes[5]--;
    nodes[6]--;
    nodes[7]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "new hex");
    RES(cell, new_cell, "hex index");
  }

  fclose(file);

  return REF_SUCCESS;
}
static REF_STATUS ref_import_surf(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                  const char *filename) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  FILE *file;
  REF_INT nnode, ntri, nqua;
  REF_INT node, new_node;
  REF_DBL xyz[3];
  REF_INT tri, new_tri;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT qua, new_qua;
  REF_INT dummy;
  char buffer[1024];

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RES(1, fscanf(file, "%d", &ntri), "ntri");
  RES(1, fscanf(file, "%d", &nqua), "nqua");
  RES(1, fscanf(file, "%d", &nnode), "nnode");

  for (node = 0; node < nnode; node++) {
    RSS(ref_node_add(ref_node, node, &new_node), "new_node");
    RES(node, new_node, "node index");
    RES(1, fscanf(file, "%lf", &(xyz[0])), "x");
    RES(1, fscanf(file, "%lf", &(xyz[1])), "y");
    RES(1, fscanf(file, "%lf", &(xyz[2])), "z");
    ref_node_xyz(ref_node, 0, new_node) = xyz[0];
    ref_node_xyz(ref_node, 1, new_node) = xyz[1];
    ref_node_xyz(ref_node, 2, new_node) = xyz[2];
    fgets(buffer, sizeof(buffer), file);
  }

  RSS(ref_node_initialize_n_global(ref_node, nnode), "init glob");

  ref_cell = ref_grid_tri(ref_grid);
  nodes[3] = REF_EMPTY;
  for (tri = 0; tri < ntri; tri++) {
    for (node = 0; node < 4; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "tri");
    for (node = 0; node < 2; node++)
      RES(1, fscanf(file, "%d", &(dummy)), "extra 0 1");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_tri), "new tri");
    RES(tri, new_tri, "tri index");
  }

  ref_cell = ref_grid_qua(ref_grid);
  nodes[4] = REF_EMPTY;
  for (qua = 0; qua < nqua; qua++) {
    for (node = 0; node < 5; node++)
      RES(1, fscanf(file, "%d", &(nodes[node])), "qua");
    for (node = 0; node < 2; node++)
      RES(1, fscanf(file, "%d", &(dummy)), "extra 0 1");
    nodes[0]--;
    nodes[1]--;
    nodes[2]--;
    nodes[3]--;
    RSS(ref_cell_add(ref_cell, nodes, &new_qua), "new qua");
    RES(qua, new_qua, "qua index");
  }

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_import_bin_ugrid_chunk(FILE *file, REF_BOOL swap,
                                             REF_BOOL fat, REF_INT n,
                                             REF_INT *chunk) {
  REF_INT i;
  if (fat) {
    REF_LONG *actual;
    ref_malloc(actual, n, REF_LONG);
    REIS(n, fread(actual, sizeof(REF_LONG), (size_t)(n), file), "long chunk");
    if (swap) {
      for (i = 0; i < n; i++) {
        SWAP_LONG(actual[i]);
      }
    }
    for (i = 0; i < n; i++) {
      chunk[i] = (REF_INT)actual[i];
    }
    ref_free(actual);
  } else {
    REIS(n, fread(chunk, sizeof(REF_INT), (size_t)(n), file), "int chunk");
    if (swap) {
      for (i = 0; i < n; i++) {
        SWAP_INT(chunk[i]);
      }
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_import_bin_ugrid_c2n(REF_CELL ref_cell, REF_INT ncell,
                                           FILE *file, REF_BOOL swap,
                                           REF_BOOL fat) {
  REF_INT node_per, max_chunk, nread, chunk, cell, node, new_cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT *c2n;

  if (0 < ncell) {
    /* to skip boundary tag */
    if (ref_cell_size_per(ref_cell) > ref_cell_node_per(ref_cell))
      nodes[ref_cell_node_per(ref_cell)] = REF_EMPTY;
    node_per = ref_cell_node_per(ref_cell);
    max_chunk = MIN(1000000, ncell);
    ref_malloc(c2n, node_per * max_chunk, REF_INT);
    nread = 0;
    while (nread < ncell) {
      chunk = MIN(max_chunk, ncell - nread);
      RSS(ref_import_bin_ugrid_chunk(file, swap, fat, node_per * chunk, c2n),
          "c2n");
      for (cell = 0; cell < chunk; cell++) {
        for (node = 0; node < node_per; node++) {
          nodes[node] = c2n[node + node_per * cell] - 1;
        }
        RSS(ref_cell_add(ref_cell, nodes, &new_cell), "new cell");
        RES(cell + nread, new_cell, "cell index");
      }
      nread += chunk;
    }
    ref_free(c2n);
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_import_bin_ugrid_bound_tag(REF_CELL ref_cell,
                                                 REF_INT ncell, FILE *file,
                                                 REF_BOOL swap, REF_BOOL fat) {
  REF_INT node_per, max_chunk, nread, chunk, cell;
  REF_INT *tag;

  if (0 < ncell) {
    node_per = ref_cell_node_per(ref_cell);
    max_chunk = MIN(1000000, ncell);
    ref_malloc(tag, max_chunk, REF_INT);
    nread = 0;
    while (nread < ncell) {
      chunk = MIN(max_chunk, ncell - nread);
      RSS(ref_import_bin_ugrid_chunk(file, swap, fat, chunk, tag), "tag");
      for (cell = 0; cell < chunk; cell++) {
        ref_cell_c2n(ref_cell, node_per, nread + cell) = tag[cell];
      }
      nread += chunk;
    }
    ref_free(tag);
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_import_bin_ugrid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                       const char *filename, REF_BOOL swap,
                                       REF_BOOL fat) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;

  REF_INT node, new_node;
  REF_INT max_chunk, nread, chunk, ixyz;
  REF_DBL *xyz;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RSS(ref_import_bin_ugrid_chunk(file, swap, fat, 1, &nnode), "nnode");
  RSS(ref_import_bin_ugrid_chunk(file, swap, fat, 1, &ntri), "ntri");
  RSS(ref_import_bin_ugrid_chunk(file, swap, fat, 1, &nqua), "nqua");
  RSS(ref_import_bin_ugrid_chunk(file, swap, fat, 1, &ntet), "ntet");
  RSS(ref_import_bin_ugrid_chunk(file, swap, fat, 1, &npyr), "npyr");
  RSS(ref_import_bin_ugrid_chunk(file, swap, fat, 1, &npri), "npri");
  RSS(ref_import_bin_ugrid_chunk(file, swap, fat, 1, &nhex), "nhex");

  /* large block reads recommended for IO performance */
  max_chunk = MIN(1000000, nnode);
  ref_malloc(xyz, 3 * max_chunk, REF_DBL);
  nread = 0;
  while (nread < nnode) {
    chunk = MIN(max_chunk, nnode - nread);
    REIS(3 * chunk, fread(xyz, sizeof(REF_DBL), (size_t)(3 * chunk), file),
         "xyz");
    if (swap)
      for (ixyz = 0; ixyz < 3 * chunk; ixyz++) SWAP_DBL(xyz[ixyz]);
    for (node = 0; node < chunk; node++) {
      RSS(ref_node_add(ref_node, node + nread, &new_node), "new_node");
      ref_node_xyz(ref_node, 0, new_node) = xyz[0 + 3 * node];
      ref_node_xyz(ref_node, 1, new_node) = xyz[1 + 3 * node];
      ref_node_xyz(ref_node, 2, new_node) = xyz[2 + 3 * node];
    }
    nread += chunk;
  }
  ref_free(xyz);

  RSS(ref_node_initialize_n_global(ref_node, nnode), "init glob");

  RSS(ref_import_bin_ugrid_c2n(ref_grid_tri(ref_grid), ntri, file, swap, fat),
      "tri face nodes");
  RSS(ref_import_bin_ugrid_c2n(ref_grid_qua(ref_grid), nqua, file, swap, fat),
      "qua face nodes");

  RSS(ref_import_bin_ugrid_bound_tag(ref_grid_tri(ref_grid), ntri, file, swap,
                                     fat),
      "tri face tags");
  RSS(ref_import_bin_ugrid_bound_tag(ref_grid_qua(ref_grid), nqua, file, swap,
                                     fat),
      "tri face tags");

  RSS(ref_import_bin_ugrid_c2n(ref_grid_tet(ref_grid), ntet, file, swap, fat),
      "tet face nodes");
  RSS(ref_import_bin_ugrid_c2n(ref_grid_pyr(ref_grid), npyr, file, swap, fat),
      "pyr face nodes");
  RSS(ref_import_bin_ugrid_c2n(ref_grid_pri(ref_grid), npri, file, swap, fat),
      "pri face nodes");
  RSS(ref_import_bin_ugrid_c2n(ref_grid_hex(ref_grid), nhex, file, swap, fat),
      "hex face nodes");

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_import_r8_ugrid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                      const char *filename) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  FILE *file;
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_INT node, new_node;
  REF_DBL swapped_dbl;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tri, qua, new_tri, new_qua;
  REF_INT face_id;
  REF_INT node_per, cell, new_cell;
  REF_INT fortran_record_size;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RES(1, fread(&fortran_record_size, sizeof(REF_INT), 1, file), "nnode");
  SWAP_INT(fortran_record_size);
  REIS(7 * 4, fortran_record_size, "header start record size");

  RES(1, fread(&nnode, sizeof(REF_INT), 1, file), "nnode");
  RES(1, fread(&ntri, sizeof(REF_INT), 1, file), "ntri");
  RES(1, fread(&nqua, sizeof(REF_INT), 1, file), "nqua");
  RES(1, fread(&ntet, sizeof(REF_INT), 1, file), "ntet");
  RES(1, fread(&npyr, sizeof(REF_INT), 1, file), "npyr");
  RES(1, fread(&npri, sizeof(REF_INT), 1, file), "npri");
  RES(1, fread(&nhex, sizeof(REF_INT), 1, file), "nhex");

  RES(1, fread(&fortran_record_size, sizeof(REF_INT), 1, file), "nnode");
  SWAP_INT(fortran_record_size);
  REIS(7 * 4, fortran_record_size, "header end record size");

  SWAP_INT(nnode);
  SWAP_INT(ntri);
  SWAP_INT(nqua);
  SWAP_INT(ntet);
  SWAP_INT(npyr);
  SWAP_INT(npri);
  SWAP_INT(nhex);

  RES(1, fread(&fortran_record_size, sizeof(REF_INT), 1, file), "nnode");
  SWAP_INT(fortran_record_size);
  REIS(nnode * 3 * 8 + ntri * 4 * 4 + nqua * 5 * 4 + ntet * 4 * 4 +
           npyr * 5 * 4 + npri * 6 * 4 + nhex * 8 * 4,
       fortran_record_size, "block start record size");

  for (node = 0; node < nnode; node++) {
    RSS(ref_node_add(ref_node, node, &new_node), "new_node");
    RES(node, new_node, "node index");
    RES(1, fread(&swapped_dbl, sizeof(REF_DBL), 1, file), "x");
    SWAP_DBL(swapped_dbl);
    ref_node_xyz(ref_node, 0, new_node) = swapped_dbl;
    RES(1, fread(&swapped_dbl, sizeof(REF_DBL), 1, file), "y");
    SWAP_DBL(swapped_dbl);
    ref_node_xyz(ref_node, 1, new_node) = swapped_dbl;
    RES(1, fread(&swapped_dbl, sizeof(REF_DBL), 1, file), "z");
    SWAP_DBL(swapped_dbl);
    ref_node_xyz(ref_node, 2, new_node) = swapped_dbl;
  }

  RSS(ref_node_initialize_n_global(ref_node, nnode), "init glob");

  ref_cell = ref_grid_tri(ref_grid);
  nodes[3] = REF_EMPTY;
  for (tri = 0; tri < ntri; tri++) {
    node_per = ref_cell_node_per(ref_cell);
    for (node = 0; node < node_per; node++) {
      RES(1, fread(&(nodes[node]), sizeof(REF_INT), 1, file), "tri");
      SWAP_INT(nodes[node]);
      nodes[node]--;
    }
    RSS(ref_cell_add(ref_cell, nodes, &new_tri), "new tri");
    RES(tri, new_tri, "tri index");
  }

  ref_cell = ref_grid_qua(ref_grid);
  nodes[4] = REF_EMPTY;
  for (qua = 0; qua < nqua; qua++) {
    node_per = ref_cell_node_per(ref_cell);
    for (node = 0; node < node_per; node++) {
      RES(1, fread(&(nodes[node]), sizeof(REF_INT), 1, file), "qua");
      SWAP_INT(nodes[node]);
      nodes[node]--;
    }
    RSS(ref_cell_add(ref_cell, nodes, &new_qua), "new qua");
    RES(qua, new_qua, "qua index");
  }

  ref_cell = ref_grid_tri(ref_grid);
  for (tri = 0; tri < ntri; tri++) {
    RES(1, fread(&face_id, sizeof(REF_INT), 1, file), "tri");
    SWAP_INT(face_id);
    ref_cell_c2n(ref_cell, 3, tri) = face_id;
  }

  ref_cell = ref_grid_qua(ref_grid);
  for (qua = 0; qua < nqua; qua++) {
    RES(1, fread(&face_id, sizeof(REF_INT), 1, file), "qua");
    SWAP_INT(face_id);
    ref_cell_c2n(ref_cell, 4, qua) = face_id;
  }

  ref_cell = ref_grid_tet(ref_grid);
  for (cell = 0; cell < ntet; cell++) {
    node_per = ref_cell_node_per(ref_cell);
    for (node = 0; node < node_per; node++) {
      RES(1, fread(&(nodes[node]), sizeof(REF_INT), 1, file), "tet");
      SWAP_INT(nodes[node]);
      nodes[node]--;
    }
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "tet cell");
    RES(cell, new_cell, "tet index");
  }

  ref_cell = ref_grid_pyr(ref_grid);
  for (cell = 0; cell < npyr; cell++) {
    node_per = ref_cell_node_per(ref_cell);
    for (node = 0; node < node_per; node++) {
      RES(1, fread(&(nodes[node]), sizeof(REF_INT), 1, file), "pyr");
      SWAP_INT(nodes[node]);
      nodes[node]--;
    }
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "pyr cell");
    RES(cell, new_cell, "pyr index");
  }

  ref_cell = ref_grid_pri(ref_grid);
  for (cell = 0; cell < npri; cell++) {
    node_per = ref_cell_node_per(ref_cell);
    for (node = 0; node < node_per; node++) {
      RES(1, fread(&(nodes[node]), sizeof(REF_INT), 1, file), "pri");
      SWAP_INT(nodes[node]);
      nodes[node]--;
    }
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "pri cell");
    RES(cell, new_cell, "pri index");
  }

  ref_cell = ref_grid_hex(ref_grid);
  for (cell = 0; cell < nhex; cell++) {
    node_per = ref_cell_node_per(ref_cell);
    for (node = 0; node < node_per; node++) {
      RES(1, fread(&(nodes[node]), sizeof(REF_INT), 1, file), "hex");
      SWAP_INT(nodes[node]);
      nodes[node]--;
    }
    RSS(ref_cell_add(ref_cell, nodes, &new_cell), "hex cell");
    RES(cell, new_cell, "hex index");
  }

  RES(1, fread(&fortran_record_size, sizeof(REF_INT), 1, file), "nnode");
  SWAP_INT(fortran_record_size);
  REIS(nnode * 3 * 8 + ntri * 4 * 4 + nqua * 5 * 4 + ntet * 4 * 4 +
           npyr * 5 * 4 + npri * 6 * 4 + nhex * 8 * 4,
       fortran_record_size, "block end record size");

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_import_su2(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                 const char *filename) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  FILE *file;
  char line[1024];
  char *location;
  REF_INT ndime, npoin, nelem, nmark;
  REF_DBL x, y, z;
  REF_INT node, new_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER], ncell, cell, new_cell, cell_type;
  REF_INT faceid;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  ndime = -1;
  while (!feof(file)) {
    if (line != fgets(line, 1024, file)) return REF_SUCCESS;
    location = strchr(line, '=');
    if (NULL == location) continue;
    if (NULL != strstr(line, "NDIME")) {
      ndime = atoi(location + 1);
      printf("NDIME %d\n", ndime);
      if (2 == ndime) {
        ref_grid_twod(ref_grid) = REF_TRUE;
      }
    }
    if (NULL != strstr(line, "NPOIN")) {
      npoin = atoi(location + 1);
      printf("NPOIN %d\n", npoin);
      for (node = 0; node < npoin; node++) {
        RAS(line == fgets(line, 1024, file), "unable to read point xyz line");
        if (2 == ndime) {
          REIS(2, sscanf(line, "%lf %lf", &x, &y), "parse xyz");
          RSS(ref_node_add(ref_node, node, &new_node), "add node");
          ref_node_xyz(ref_node, 0, new_node) = x;
          ref_node_xyz(ref_node, 1, new_node) = y;
          ref_node_xyz(ref_node, 2, new_node) = 0.0;
        } else {
          REIS(3, sscanf(line, "%lf %lf %lf", &x, &y, &z), "parse xyz");
          RSS(ref_node_add(ref_node, node, &new_node), "add node");
          ref_node_xyz(ref_node, 0, new_node) = x;
          ref_node_xyz(ref_node, 1, new_node) = y;
          ref_node_xyz(ref_node, 2, new_node) = z;
        }
      }
      RSS(ref_node_initialize_n_global(ref_node, npoin), "init glob");
    }
    if (NULL != strstr(line, "NELEM")) {
      nelem = atoi(location + 1);
      printf("NELEM %d\n", nelem);
      for (cell = 0; cell < nelem; cell++) {
        RAS(line == fgets(line, 1024, file), "unable to read element line");
        REIS(1, sscanf(line, "%d", &cell_type), "parse element type");
        switch (cell_type) {
          case VTK_TETRA:
            REIS(5,
                 sscanf(line, "%d %d %d %d %d", &cell_type, &(nodes[0]),
                        &(nodes[1]), &(nodes[2]), &(nodes[3])),
                 "parse element");
            RSS(ref_cell_add(ref_grid_tet(ref_grid), nodes, &new_cell), "tet");
            break;
          case VTK_PYRAMID:
            REIS(6,
                 sscanf(line, "%d %d %d %d %d %d", &cell_type, &(nodes[0]),
                        &(nodes[1]), &(nodes[2]), &(nodes[3]), &(nodes[4])),
                 "parse element");
            VTK_PYRAMID_TO_UGRID(nodes);
            RSS(ref_cell_add(ref_grid_pyr(ref_grid), nodes, &new_cell), "pyr");
            break;
          case VTK_WEDGE:
            REIS(7,
                 sscanf(line, "%d %d %d %d %d %d %d", &cell_type, &(nodes[0]),
                        &(nodes[1]), &(nodes[2]), &(nodes[3]), &(nodes[4]),
                        &(nodes[5])),
                 "parse element");
            VTK_WEDGE_TO_UGRID(nodes);
            RSS(ref_cell_add(ref_grid_pri(ref_grid), nodes, &new_cell), "tri");
            break;
          case VTK_HEXAHEDRON:
            REIS(9,
                 sscanf(line, "%d %d %d %d %d %d %d %d %d", &cell_type,
                        &(nodes[0]), &(nodes[1]), &(nodes[2]), &(nodes[3]),
                        &(nodes[4]), &(nodes[5]), &(nodes[6]), &(nodes[7])),
                 "parse element");
            RSS(ref_cell_add(ref_grid_hex(ref_grid), nodes, &new_cell), "tri");
            break;
          case VTK_TRIANGLE:
            REIS(4,
                 sscanf(line, "%d %d %d %d", &cell_type, &(nodes[0]),
                        &(nodes[1]), &(nodes[2])),
                 "parse marker element");
            nodes[3] = 0;
            RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &new_cell), "tri");
            break;
          case VTK_QUAD:
            REIS(5,
                 sscanf(line, "%d %d %d %d %d", &cell_type, &(nodes[0]),
                        &(nodes[1]), &(nodes[2]), &(nodes[3])),
                 "parse marker element");
            nodes[4] = 0;
            RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &new_cell), "tri");
            break;

          default:
            printf("cell_type = %d\n", cell_type);
            THROW("unknown SU2/VTK ELEM type");
        }
      }
    }
    if (NULL != strstr(line, "NMARK")) {
      nmark = atoi(location + 1);
      printf("NMARK %d\n", nmark);
      for (faceid = 1; faceid <= nmark; faceid++) {
        RAS(line == fgets(line, 1024, file), "unable to read marker tag");
        printf("%d: %s", faceid, line);
        RAS(NULL != strstr(line, "MARKER_TAG"), "MARKER_TAG not found");
        RAS(line == fgets(line, 1024, file), "unable to read marker elems");
        RAS(NULL != strstr(line, "MARKER_ELEMS"), "MARKER_ELEMS not found");
        location = strchr(line, '=');
        RNS(location, "MARKER_ELEMS missing =");
        ncell = atoi(location + 1);
        printf("%d: %s", ncell, line);
        for (cell = 0; cell < ncell; cell++) {
          RAS(line == fgets(line, 1024, file), "unable to read marker line");
          REIS(1, sscanf(line, "%d", &cell_type), "parse marker element type");
          switch (cell_type) {
            case VTK_TRIANGLE:
              REIS(4,
                   sscanf(line, "%d %d %d %d", &cell_type, &(nodes[0]),
                          &(nodes[1]), &(nodes[2])),
                   "parse marker element");
              nodes[3] = faceid;
              RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &new_cell),
                  "tri");
              break;
            case VTK_QUAD:
              REIS(5,
                   sscanf(line, "%d %d %d %d %d", &cell_type, &(nodes[0]),
                          &(nodes[1]), &(nodes[2]), &(nodes[3])),
                   "parse marker element");
              nodes[4] = faceid;
              RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &new_cell),
                  "tri");
              break;
            case VTK_LINE:
              REIS(3,
                   sscanf(line, "%d %d %d", &cell_type, &(nodes[0]),
                          &(nodes[1])),
                   "parse marker element");
              nodes[2] = faceid;
              RSS(ref_cell_add(ref_grid_edg(ref_grid), nodes, &new_cell),
                  "edg");
              break;
            default:
              printf("cell_type = %d\n", cell_type);
              THROW("unknown SU2/VTK MARKER ELEM type");
          }
        }
      }
    }
  }

  fclose(file);

  return REF_SUCCESS;
}

/*
                           2
                        ./   \
                     .  / nt2 \
                  .    /       \
               .      /nt0  nt1 \
            .        0-ne0---ne1-1  y=1 second plane, face 1
         5       .           .
       /   \  .           .
      /  t1 \          .
     /   .   \      .
    / t0   t2 \  .
   3--e0----e1-4   y=0 first plane, face 2

 */

static REF_STATUS ref_import_msh(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                 const char *filename) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  FILE *file;
  char line[1024];
  REF_INT dummy, row;
  REF_DBL x, y, z;
  REF_INT dim, nnode, node, new_node;
  REF_INT nedge, edge, n0, n1, n2, n3, id;
  REF_INT ntri, tri;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER], new_cell;
  REF_INT status;
  REF_INT elem, nelem, type, flag, three, zero;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  while (!feof(file)) {
    status = fscanf(file, "%s", line);
    if (EOF == status) return REF_SUCCESS;
    REIS(1, status, "line read failed");

    if (0 == strcmp("Dimension", line)) {
      REIS(1, fscanf(file, "%d", &dim), "read dim");
      if (2 == dim) ref_grid_twod(ref_grid) = REF_TRUE;
    }

    if (0 == strcmp("Vertices", line)) {
      REIS(1, fscanf(file, "%d", &nnode), "read nnode");
      for (node = 0; node < nnode; node++) {
        REIS(3, fscanf(file, "%lf %lf %d", &x, &y, &dummy), "read xy");
        RSS(ref_node_add(ref_node, node, &new_node), "add node");
        ref_node_xyz(ref_node, 0, new_node) = x;
        ref_node_xyz(ref_node, 1, new_node) = y;
        ref_node_xyz(ref_node, 2, new_node) = 0.0;
      }
      RSS(ref_node_initialize_n_global(ref_node, nnode), "init glob");
    }

    if (0 == strcmp("Edges", line)) {
      REIS(1, fscanf(file, "%d", &nedge), "read nedge");
      for (edge = 0; edge < nedge; edge++) {
        REIS(3, fscanf(file, "%d %d %d", &n0, &n1, &id), "read edge");
        n0--;
        n1--;
        nodes[0] = n0;
        nodes[1] = n1;
        nodes[2] = id;
        RSS(ref_cell_add(ref_grid_edg(ref_grid), nodes, &new_cell),
            "quad face for an edge");
      }
    }

    if (0 == strcmp("Triangles", line)) {
      REIS(1, fscanf(file, "%d", &ntri), "read ntri");
      for (tri = 0; tri < ntri; tri++) {
        REIS(4, fscanf(file, "%d %d %d %d", &n0, &n1, &n2, &id), "read tri");
        n0--;
        n1--;
        n2--;
        nodes[0] = n0;
        nodes[1] = n1;
        nodes[2] = n2;
        nodes[3] = id;
        RSS(ref_cell_add(ref_grid_tri(ref_grid), nodes, &new_cell),
            "tri face for tri");
      }
    }

    if (0 == strcmp("Quadrilaterals", line)) {
      REIS(1, fscanf(file, "%d", &ntri), "read ntri");
      for (tri = 0; tri < ntri; tri++) {
        REIS(5, fscanf(file, "%d %d %d %d %d", &n0, &n1, &n2, &n3, &id),
             "read quad");
        n0--;
        n1--;
        n2--;
        n3--;
        nodes[0] = n0;
        nodes[1] = n1;
        nodes[2] = n2;
        nodes[3] = n3;
        nodes[4] = id;
        RSS(ref_cell_add(ref_grid_qua(ref_grid), nodes, &new_cell),
            "qua face for qua");
      }
    }

    if (0 == strcmp("$Nodes", line)) {
      REIS(1, fscanf(file, "%d", &nnode), "read nnode");
      printf("$Nodes\n%d\n", nnode);
      for (node = 0; node < nnode; node++) {
        REIS(4, fscanf(file, "%d %lf %lf %lf", &row, &x, &y, &z),
             "read $Nodes xyz");
        REIS(node + 1, row, "row index miss match in $Nodes");
        RSS(ref_node_add(ref_node, node, &new_node), "add node");
        ref_node_xyz(ref_node, 0, new_node) = x;
        ref_node_xyz(ref_node, 1, new_node) = y;
        ref_node_xyz(ref_node, 2, new_node) = z;
      }
      RSS(ref_node_initialize_n_global(ref_node, nnode), "init glob");
    }

    if (0 == strcmp("$Elements", line)) {
      REIS(1, fscanf(file, "%d", &nelem), "read nelements");
      printf("$Elements\n%d\n", nelem);
      for (elem = 0; elem < nelem; elem++) {
        REIS(6,
             fscanf(file, "%d %d %d %d %d %d", &row, &type, &three, &id, &flag,
                    &zero),
             "$Elements description");
        REIS(elem + 1, row, "row index miss match in $Elements");
        switch (type) {
          case 5:
            ref_cell = ref_grid_hex(ref_grid);
            break;
          case 3:
            ref_cell = ref_grid_qua(ref_grid);
            break;
          default:
            printf("type = %d\n", type);
            THROW("unknown $Elements type");
        }
        for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
          REIS(1, fscanf(file, "%d", &(nodes[node])), "$Elements node");
          (nodes[node])--;
        }
        if (3 == type) nodes[4] = id;
        RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add $Element");
      }
      ref_grid_inspect(ref_grid);
    }
  }

  fclose(file);

  return REF_IMPLEMENT;
}

static REF_STATUS meshb_real(FILE *file, REF_INT version, REF_DBL *real) {
  float temp_float;
  double temp_double;

  if (1 == version) {
    REIS(1, fread(&temp_float, sizeof(temp_float), 1, file), "read float");
    *real = (REF_DBL)temp_float;
  } else {
    REIS(1, fread(&temp_double, sizeof(temp_double), 1, file), "read double");
    *real = (REF_DBL)temp_double;
  }

  return REF_SUCCESS;
}

static REF_STATUS meshb_pos(FILE *file, REF_INT version, REF_FILEPOS *pos) {
  int temp_int;
  long temp_long;

  if (3 <= version) {
    REIS(1, fread(&temp_long, sizeof(temp_long), 1, file), "read long");
    *pos = temp_long;
  } else {
    REIS(1, fread(&temp_int, sizeof(temp_int), 1, file), "read double");
    *pos = temp_int;
  }

  return REF_SUCCESS;
}
REF_STATUS ref_import_meshb_header(const char *filename, REF_INT *version,
                                   REF_FILEPOS *key_pos) {
  FILE *file;
  int int_code, int_version;
  REF_INT keyword_code;
  REF_FILEPOS position, next_position, end_position;

  for (keyword_code = 0; keyword_code < REF_IMPORT_MESHB_LAST_KEYWORD;
       keyword_code++)
    key_pos[keyword_code] = REF_EMPTY;

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  REIS(1, fread((unsigned char *)&int_code, 4, 1, file), "code");
  REIS(1, int_code, "code");
  REIS(1, fread((unsigned char *)&int_version, 4, 1, file), "version");
  if (int_version < 1 || 4 < int_version) {
    printf("version %d not supported\n", int_version);
    THROW("version");
  }
  *version = (REF_INT)int_version;

  next_position = ftello(file);
  REIS(0, fseeko(file, 0, SEEK_END), "fseeko END failed");
  end_position = ftello(file);
  while (next_position <= end_position && 0 != next_position) {
    position = next_position;
    REIS(0, fseeko(file, position, SEEK_SET), "fseeko NEXT failed");
    /*
#include <inttypes.h>
    , {
      printf("end_position = %jd\n", (intmax_t)end_position);
      printf("next_position = %jd\n", (intmax_t)next_position);
      for (keyword_code = 0; keyword_code < REF_IMPORT_MESHB_LAST_KEYWORD;
           keyword_code++) {
        if (REF_EMPTY != key_pos[keyword_code]) {
          printf("key_pos[%d] = %jd\n", keyword_code,
                 (intmax_t)key_pos[keyword_code]);
        }
      }
    });
    */
    REIS(1, fread((unsigned char *)&keyword_code, 4, 1, file), "keyword code");
    if (0 <= keyword_code && keyword_code < REF_IMPORT_MESHB_LAST_KEYWORD) {
      key_pos[keyword_code] = position;
    } else {
      printf("ignoring keyword %d\n", keyword_code);
    }
    RSS(meshb_pos(file, *version, &next_position), "pos");
  }

  fclose(file);
  return REF_SUCCESS;
}

REF_STATUS ref_import_meshb_jump(FILE *file, REF_INT version,
                                 REF_FILEPOS *key_pos, REF_INT keyword,
                                 REF_BOOL *available,
                                 REF_FILEPOS *next_position) {
  REF_INT keyword_code;
  REF_FILEPOS position;

  if (keyword < 0 && REF_IMPORT_MESHB_LAST_KEYWORD <= keyword) {
    *available = REF_FALSE;
    *next_position = 0;
    return REF_INVALID;
  }
  position = key_pos[keyword];

  if ((REF_FILEPOS)REF_EMPTY != position) {
    *available = REF_TRUE;
  } else {
    *available = REF_FALSE;
    *next_position = 0;
    return REF_SUCCESS;
  }
  REIS(0, fseeko(file, position, SEEK_SET), "fseeko keyword failed");
  REIS(1, fread((unsigned char *)&keyword_code, 4, 1, file), "keyword code");
  REIS(keyword, keyword_code, "keyword code");
  RSS(meshb_pos(file, version, next_position), "pos");
  return REF_SUCCESS;
}

static REF_STATUS ref_import_meshb_int(FILE *file, REF_INT version,
                                       REF_INT *value) {
  int int_value;
  long long_value;
  if (version < 4) {
    REIS(1, fread(&int_value, sizeof(int), 1, file), "int value");
    *value = int_value;
  } else {
    REIS(1, fread(&long_value, sizeof(long), 1, file), "long value");
    *value = (REF_INT)long_value;
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_import_meshb_size(FILE *file, REF_INT version,
                                        REF_SIZE *value) {
  unsigned int int_value;
  unsigned long long_value;
  if (version < 4) {
    REIS(1, fread(&int_value, sizeof(int), 1, file), "int value");
    *value = (REF_SIZE)int_value;
  } else {
    REIS(1, fread(&long_value, sizeof(long), 1, file), "long value");
    *value = (REF_SIZE)long_value;
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_import_meshb(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                   const char *filename) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_GEOM ref_geom;
  FILE *file;
  REF_INT version, dim;
  REF_BOOL available;
  REF_FILEPOS next_position;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  REF_INT keyword_code, nnode, node, new_node;
  REF_INT ncell, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER + 1]; /* everyone gets id in meshb */
  REF_INT new_cell;
  REF_INT n0, n1, n2, n3, n4, id, group, node_per;
  REF_INT geom_keyword, type, i, geom, ngeom;
  REF_DBL param[2];
  REF_INT cad_data_keyword;
  REF_BOOL verbose = REF_FALSE;

  if (verbose) printf("header %s\n", filename);
  RSS(ref_import_meshb_header(filename, &version, key_pos), "header");
  if (verbose) printf("meshb version %d\n", version);

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);
  ref_geom = ref_grid_geom(ref_grid);

  if (verbose) printf("open %s\n", filename);
  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  RSS(ref_import_meshb_jump(file, version, key_pos, 3, &available,
                            &next_position),
      "jump");
  RAS(available, "meshb missing dimension");
  REIS(1, fread((unsigned char *)&dim, 4, 1, file), "dim");
  if (verbose) printf("meshb dim %d\n", dim);
  if (dim < 2 || 3 < dim) {
    printf("dim %d not supported\n", dim);
    THROW("dim");
  }
  if (2 == dim) ref_grid_twod(ref_grid) = REF_TRUE;

  RSS(ref_import_meshb_jump(file, version, key_pos, 4, &available,
                            &next_position),
      "jump");
  RAS(available, "meshb missing vertex");
  RSS(ref_import_meshb_int(file, version, &nnode), "nnode");
  if (verbose) printf("nnode %d\n", nnode);

  for (node = 0; node < nnode; node++) {
    RSS(ref_node_add(ref_node, node, &new_node), "add node");
    if (2 == dim) {
      RSS(meshb_real(file, version, &(ref_node_xyz(ref_node, 0, new_node))),
          "x");
      RSS(meshb_real(file, version, &(ref_node_xyz(ref_node, 1, new_node))),
          "y");
      ref_node_xyz(ref_node, 2, new_node) = 0.0;
    } else {
      RSS(meshb_real(file, version, &(ref_node_xyz(ref_node, 0, new_node))),
          "x");
      RSS(meshb_real(file, version, &(ref_node_xyz(ref_node, 1, new_node))),
          "y");
      RSS(meshb_real(file, version, &(ref_node_xyz(ref_node, 2, new_node))),
          "z");
    }
    RSS(ref_import_meshb_int(file, version, &id), "nnode");
  }
  REIS(next_position, ftello(file), "end location");

  RSS(ref_node_initialize_n_global(ref_node, nnode), "init glob");

  each_ref_grid_all_ref_cell(ref_grid, group, ref_cell) {
    RSS(ref_cell_meshb_keyword(ref_cell, &keyword_code), "kw");
    RSS(ref_import_meshb_jump(file, version, key_pos, keyword_code, &available,
                              &next_position),
        "jump");
    if (available) {
      node_per = ref_cell_node_per(ref_cell);
      RSS(ref_import_meshb_int(file, version, &ncell), "ncell");
      if (verbose) printf(" group %d ncell %d\n", group, ncell);
      for (cell = 0; cell < ncell; cell++) {
        for (node = 0; node < (1 + node_per); node++) {
          RSS(ref_import_meshb_int(file, version, &(nodes[node])), "c2n");
        }
        for (node = 0; node < node_per; node++) {
          nodes[node]--;
        }
        if (REF_CELL_PYR == ref_cell_type(ref_cell)) {
          /* convention: square basis is 0-1-2-3
             (oriented counter clockwise like trias) and top vertex is 4 */
          n0 = nodes[0];
          n1 = nodes[1];
          n2 = nodes[2];
          n3 = nodes[3];
          n4 = nodes[4];
          nodes[0] = n0;
          nodes[3] = n1;
          nodes[4] = n2;
          nodes[1] = n3;
          nodes[2] = n4;
        }
        RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add cell");
      }
      REIS(next_position, ftell(file), "cell inconsistent");
    }
  }

  each_ref_type(ref_geom, type) {
    geom_keyword = 40 + type;
    RSS(ref_import_meshb_jump(file, version, key_pos, geom_keyword, &available,
                              &next_position),
        "jump");
    if (available) {
      RSS(ref_import_meshb_int(file, version, &(ngeom)), "ngeom");
      if (verbose) printf("type %d ngeom %d\n", type, ngeom);

      for (geom = 0; geom < ngeom; geom++) {
        RSS(ref_import_meshb_int(file, version, &(node)), "node");
        RSS(ref_import_meshb_int(file, version, &(id)), "node");
        for (i = 0; i < type; i++)
          REIS(1, fread(&(param[i]), sizeof(double), 1, file), "param");
        node--;
        RSS(ref_geom_add(ref_geom, node, type, id, param), "add geom");
        if (0 < type) {
          double double_gref;
          REF_INT new_geom;
          REIS(1, fread(&(double_gref), sizeof(double), 1, file), "gref");
          RSS(ref_geom_find(ref_geom, node, type, id, &new_geom), "find");
          ref_geom_gref(ref_geom, new_geom) = (REF_INT)double_gref;
        }
      }
      REIS(next_position, ftello(file), "end location");
    }
  }

  cad_data_keyword = 126; /* GmfByteFlow */
  RSS(ref_import_meshb_jump(file, version, key_pos, cad_data_keyword,
                            &available, &next_position),
      "jump");
  if (available) {
    RSS(ref_import_meshb_size(file, version,
                              &(ref_geom_cad_data_size(ref_geom))),
        "cad data size");
    if (verbose)
      printf("cad_data %ld bytes\n", (long)ref_geom_cad_data_size(ref_geom));
    /* safe non-NULL free, if already allocated, to prevent memory leaks */
    ref_free(ref_geom_cad_data(ref_geom));
    ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                      ref_geom_cad_data_size(ref_geom), REF_BYTE);
    REIS(ref_geom_cad_data_size(ref_geom),
         fread(ref_geom_cad_data(ref_geom), sizeof(REF_BYTE),
               ref_geom_cad_data_size(ref_geom), file),
         "cad_data");
    REIS(next_position, ftello(file), "end location");
  }

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_import_by_extension(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                   const char *filename) {
  size_t end_of_string;

  end_of_string = strlen(filename);

  if (strcmp(&filename[end_of_string - 10], ".lb8.ugrid") == 0) {
    RSS(ref_import_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_FALSE,
                             REF_FALSE),
        "lb8_ugrid failed");
  } else if (strcmp(&filename[end_of_string - 9], ".b8.ugrid") == 0) {
    RSS(ref_import_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_TRUE,
                             REF_FALSE),
        "b8_ugrid failed");
  } else if (strcmp(&filename[end_of_string - 11], ".lb8l.ugrid") == 0) {
    RSS(ref_import_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_FALSE,
                             REF_TRUE),
        "lb8_ugrid failed");
  } else if (strcmp(&filename[end_of_string - 10], ".b8l.ugrid") == 0) {
    RSS(ref_import_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_TRUE,
                             REF_TRUE),
        "b8_ugrid failed");
  } else if (strcmp(&filename[end_of_string - 12], ".lb8.ugrid64") == 0) {
    RSS(ref_import_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_FALSE,
                             REF_TRUE),
        "lb8_ugrid failed");
  } else if (strcmp(&filename[end_of_string - 11], ".b8.ugrid64") == 0) {
    RSS(ref_import_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_TRUE,
                             REF_TRUE),
        "b8_ugrid failed");
  } else if (strcmp(&filename[end_of_string - 9], ".r8.ugrid") == 0) {
    RSS(ref_import_r8_ugrid(ref_grid_ptr, ref_mpi, filename),
        "r8_ugrid failed");
  } else if (strcmp(&filename[end_of_string - 6], ".ugrid") == 0) {
    RSS(ref_import_ugrid(ref_grid_ptr, ref_mpi, filename), "ugrid failed");
  } else if (strcmp(&filename[end_of_string - 5], ".surf") == 0) {
    RSS(ref_import_surf(ref_grid_ptr, ref_mpi, filename), "surf failed");
  } else if (strcmp(&filename[end_of_string - 6], ".fgrid") == 0) {
    RSS(ref_import_fgrid(ref_grid_ptr, ref_mpi, filename), "fgrid failed");
  } else if (strcmp(&filename[end_of_string - 4], ".su2") == 0) {
    RSS(ref_import_su2(ref_grid_ptr, ref_mpi, filename), "su2 failed");
  } else if (strcmp(&filename[end_of_string - 4], ".msh") == 0) {
    RSS(ref_import_msh(ref_grid_ptr, ref_mpi, filename), "msh failed");
  } else if (strcmp(&filename[end_of_string - 6], ".meshb") == 0) {
    RSS(ref_import_meshb(ref_grid_ptr, ref_mpi, filename), "meshb failed");
  } else {
    printf("%s: %d: %s %s\n", __FILE__, __LINE__,
           "input file name extension unknown", filename);
    RSS(REF_FAILURE, "unknown file extension");
  }

  RSS(ref_grid_inward_boundary_orientation(*ref_grid_ptr),
      "inward boundary orientation");

  return REF_SUCCESS;
}

REF_STATUS ref_import_examine_header(const char *filename) {
  FILE *file;
  REF_FILEPOS next_position, end_position;
  int i4, i4_swapped, version, keyword_code;
  long i8, i8_swapped;
  int i;
  REF_BOOL file_position_report = REF_FALSE;

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  printf(" -- 32bit ugrid header\n");

  rewind(file);

  for (i = 0; i < 8; i++) {
    RES(1, fread(&i4, sizeof(int), 1, file), "int");
    i4_swapped = i4;
    SWAP_INT(i4_swapped);
    printf(" %d: %d (%d swapped) ints\n", i, i4, i4_swapped);
  }

  printf(" -- 64bit ugrid header\n");

  rewind(file);

  for (i = 0; i < 3; i++) {
    RES(1, fread(&i4, sizeof(int), 1, file), "int");
  }
  for (i = 3; i < 7; i++) {
    RES(1, fread(&i8, sizeof(long), 1, file), "long");
    i8_swapped = i8;
    SWAP_INT(i8_swapped);
    printf(" %d: %ld (%ld swapped) long\n", i, i8, i8_swapped);
  }

  printf(" -- meshb/solb\n");

  rewind(file);

  printf("%d sizeof(REF_FILEPOS)\n", (REF_INT)sizeof(REF_FILEPOS));

  REIS(0, fseeko(file, 0, SEEK_END), "fseeko END failed");
  end_position = ftello(file);
  if (file_position_report) printf("%ld end_position\n", (long)end_position);
  rewind(file);

  REIS(1, fread((unsigned char *)&i4, 4, 1, file), "code");
  printf("%d meshb code\n", i4);
  if (1 != i4) goto close_file_and_return;
  REIS(1, fread((unsigned char *)&i4, 4, 1, file), "version");
  printf("%d version\n", i4);
  if (1 > i4 || i4 > 4) goto close_file_and_return;
  version = i4;
  next_position = ftello(file);
  while (next_position <= end_position && 0 < next_position) {
    REIS(0, fseeko(file, next_position, SEEK_SET), "fseeko NEXT failed");
    if (file_position_report)
      printf("%ld current position\n", (long)next_position);
    REIS(1, fread((unsigned char *)&keyword_code, 4, 1, file), "keyword code");
    printf("%d keyword", keyword_code);
    switch (keyword_code) {
      case 3:
        printf(" dimension\n");
        break;
      case 4:
        printf(" vertex\n");
        break;
      case 5:
        printf(" edge\n");
        break;
      case 6:
        printf(" triangle\n");
        break;
      case 7:
        printf(" quad\n");
        break;
      case 8:
        printf(" tetrahedra\n");
        break;
      case 9:
        printf(" prism\n");
        break;
      case 10:
        printf(" hex\n");
        break;
      case 49:
        printf(" pyramid\n");
        break;
      case 40:
        printf(" geom node assoc\n");
        break;
      case 41:
        printf(" geom edge assoc\n");
        break;
      case 42:
        printf(" geom face assoc\n");
        break;
      case 62:
        printf(" solution at vertices\n");
        break;
      case 126:
        printf(" CAD data\n");
        break;
      case 54:
        printf(" END\n");
        break;
      default:
        printf("\n");
    }
    RSS(meshb_pos(file, version, &next_position), "meshb pos");
    if (file_position_report) printf("%ld next position", (long)next_position);
    if (next_position > 0) {
      printf(" %ld size\n", (long)(next_position - ftello(file)));
    } else {
      printf("\n");
    }
    if (version >= 4 && keyword_code != 3) {
      if (ftello(file) < end_position) {
        REIS(1, fread((unsigned char *)&i8, 8, 1, file), "code");
        printf("  %ld first i8\n", i8);
      }
    } else {
      if (ftello(file) < end_position) {
        REIS(1, fread((unsigned char *)&i4, 4, 1, file), "code");
        printf("  %d first i4\n", i4);
      }
    }
    if (ftello(file) < end_position) {
      REIS(1, fread((unsigned char *)&i4, 4, 1, file), "code");
      printf("  %d second i4\n", i4);
    }
  }

close_file_and_return:
  fclose(file);
  return REF_SUCCESS;
}
