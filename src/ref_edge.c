
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

#include <stdio.h>
#include <stdlib.h>

#include "ref_edge.h"

#include "ref_malloc.h"
#include "ref_mpi.h"

static REF_STATUS ref_edge_builder_uniq(REF_EDGE ref_edge, REF_GRID ref_grid) {
  REF_INT group, cell, cell_edge;
  REF_INT node0, node1;
  REF_CELL ref_cell;

  /* first allocation to an estimated size */
  if (NULL == (void *)ref_edge->e2n) {
    REF_INT edge_per_node_estimate = 8; /* 7 is the tet estimate */
    REIS(0, ref_edge_max(ref_edge), "should be zero size");
    ref_edge_max(ref_edge) =
        MAX(100, edge_per_node_estimate * ref_node_n(ref_grid_node(ref_grid)));
    ref_malloc_init(ref_edge->e2n, 2 * ref_edge_max(ref_edge), REF_INT,
                    REF_EMPTY);
  }

  each_ref_grid_ref_cell(ref_grid, group, ref_cell) {
    each_ref_cell_valid_cell(ref_cell, cell) {
      each_ref_cell_cell_edge(ref_cell, cell_edge) {
        node0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
        node1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
        RSS(ref_edge_uniq(ref_edge, node0, node1), "add uniq");
      }
    }
  }

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    each_ref_cell_cell_edge(ref_cell, cell_edge) {
      node0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
      node1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
      RSS(ref_edge_uniq(ref_edge, node0, node1), "add uniq");
    }
  }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell(ref_cell, cell) {
    each_ref_cell_cell_edge(ref_cell, cell_edge) {
      node0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
      node1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
      RSS(ref_edge_uniq(ref_edge, node0, node1), "add uniq");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_edge_create(REF_EDGE *ref_edge_ptr, REF_GRID ref_grid) {
  REF_EDGE ref_edge;

  ref_malloc(*ref_edge_ptr, 1, REF_EDGE_STRUCT);

  ref_edge = *ref_edge_ptr;

  ref_edge_n(ref_edge) = 0;
  ref_edge_max(ref_edge) = 0;
  ref_edge->e2n = (REF_INT *)NULL;

  RSS(ref_adj_create(&(ref_edge_adj(ref_edge))), "create adj");

  ref_edge_node(ref_edge) = ref_grid_node(ref_grid);

  RSS(ref_edge_builder_uniq(ref_edge, ref_grid), "build edges");

  return REF_SUCCESS;
}

REF_STATUS ref_edge_free(REF_EDGE ref_edge) {
  if (NULL == (void *)ref_edge) return REF_NULL;

  RSS(ref_adj_free(ref_edge_adj(ref_edge)), "free adj");
  ref_free(ref_edge->e2n);

  ref_free(ref_edge);

  return REF_SUCCESS;
}

REF_STATUS ref_edge_uniq(REF_EDGE ref_edge, REF_INT node0, REF_INT node1) {
  REF_INT edge;

  /* do nothing if we already have it */
  RXS(ref_edge_with(ref_edge, node0, node1, &edge), REF_NOT_FOUND,
      "find existing");
  if (REF_EMPTY != edge) return REF_SUCCESS;

  /* incemental reallocation */
  if (ref_edge_n(ref_edge) >= ref_edge_max(ref_edge)) {
    REF_INT orig, chunk;
    orig = ref_edge_max(ref_edge);
    /* geometric growth for efficiency */
    chunk = MAX(5000, (REF_INT)(1.5 * (REF_DBL)orig));
    ref_edge_max(ref_edge) = orig + chunk;

    ref_realloc(ref_edge->e2n, 2 * ref_edge_max(ref_edge), REF_INT);
    for (edge = orig; edge < ref_edge_max(ref_edge); edge++) {
      ref_edge_e2n(ref_edge, 0, edge) = REF_EMPTY;
      ref_edge_e2n(ref_edge, 1, edge) = REF_EMPTY;
    }
  }

  edge = ref_edge_n(ref_edge);
  ref_edge_n(ref_edge)++;
  ref_edge_e2n(ref_edge, 0, edge) = node0;
  ref_edge_e2n(ref_edge, 1, edge) = node1;

  RSS(ref_adj_add(ref_edge_adj(ref_edge), ref_edge_e2n(ref_edge, 0, edge),
                  edge),
      "adj n0");
  RSS(ref_adj_add(ref_edge_adj(ref_edge), ref_edge_e2n(ref_edge, 1, edge),
                  edge),
      "adj n1");

  return REF_SUCCESS;
}

REF_STATUS ref_edge_with(REF_EDGE ref_edge, REF_INT node0, REF_INT node1,
                         REF_INT *edge) {
  REF_INT item, ref;
  REF_INT n0, n1;

  *edge = REF_EMPTY;

  each_ref_adj_node_item_with_ref(ref_edge_adj(ref_edge), node0, item, ref) {
    n0 = ref_edge_e2n(ref_edge, 0, ref);
    n1 = ref_edge_e2n(ref_edge, 1, ref);
    if ((n0 == node0 && n1 == node1) || (n0 == node1 && n1 == node0)) {
      *edge = ref;
      return REF_SUCCESS;
    }
  }

  return REF_NOT_FOUND;
}

REF_STATUS ref_edge_part(REF_EDGE ref_edge, REF_INT edge, REF_INT *part) {
  REF_NODE ref_node = ref_edge_node(ref_edge);

  if (ref_node_global(ref_node, ref_edge_e2n(ref_edge, 0, edge)) <
      ref_node_global(ref_node, ref_edge_e2n(ref_edge, 1, edge))) {
    *part = ref_node_part(ref_node, ref_edge_e2n(ref_edge, 0, edge));
  } else {
    *part = ref_node_part(ref_node, ref_edge_e2n(ref_edge, 1, edge));
  }

  return REF_SUCCESS;
}

REF_STATUS ref_edge_ghost_int(REF_EDGE ref_edge, REF_MPI ref_mpi,
                              REF_INT *data) {
  REF_NODE ref_node = ref_edge_node(ref_edge);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT edge;
  REF_INT part;
  REF_INT *a_next, *a_edge;
  REF_GLOB *a_nodes, *b_nodes;
  REF_INT *a_data, *b_data;

  REF_INT node0, node1;
  REF_INT request;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (part != ref_mpi_rank(ref_mpi)) a_size[part]++;
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_total += a_size[part];
  ref_malloc(a_nodes, 2 * a_total, REF_GLOB);
  ref_malloc(a_data, a_total, REF_INT);
  ref_malloc(a_edge, a_total, REF_INT);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_total += b_size[part];
  ref_malloc(b_nodes, 2 * b_total, REF_GLOB);
  ref_malloc(b_data, b_total, REF_INT);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) a_next[part] =
      a_next[part - 1] + a_size[part - 1];

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (part != ref_mpi_rank(ref_mpi)) {
      a_edge[a_next[part]] = edge;
      a_nodes[0 + 2 * a_next[part]] =
          ref_node_global(ref_node, ref_edge_e2n(ref_edge, 0, edge));
      a_nodes[1 + 2 * a_next[part]] =
          ref_node_global(ref_node, ref_edge_e2n(ref_edge, 1, edge));
      (a_next[part])++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_nodes, a_size, b_nodes, b_size, 2,
                        REF_GLOB_TYPE),
      "alltoallv requested nodes");

  for (request = 0; request < b_total; request++) {
    RSS(ref_node_local(ref_node, b_nodes[0 + 2 * request], &node0), "loc 0");
    RSS(ref_node_local(ref_node, b_nodes[1 + 2 * request], &node1), "loc 1");
    RSS(ref_edge_with(ref_edge, node0, node1, &edge), "find edge");
    b_data[request] = data[edge];
  }

  RSS(ref_mpi_alltoallv(ref_mpi, b_data, b_size, a_data, a_size, 1,
                        REF_INT_TYPE),
      "alltoallv return data");

  for (request = 0; request < a_total; request++) {
    data[a_edge[request]] = a_data[request];
  }

  ref_free(a_next);

  ref_free(b_data);
  ref_free(b_nodes);

  ref_free(a_edge);

  ref_free(a_data);
  ref_free(a_nodes);

  ref_free(b_size);
  ref_free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_edge_ghost_glob(REF_EDGE ref_edge, REF_MPI ref_mpi,
                               REF_GLOB *data) {
  REF_NODE ref_node = ref_edge_node(ref_edge);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT edge;
  REF_INT part;
  REF_INT *a_next, *a_edge;
  REF_GLOB *a_nodes, *b_nodes;
  REF_GLOB *a_data, *b_data;

  REF_INT node0, node1;
  REF_INT request;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (part != ref_mpi_rank(ref_mpi)) a_size[part]++;
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_total += a_size[part];
  ref_malloc(a_nodes, 2 * a_total, REF_GLOB);
  ref_malloc(a_data, a_total, REF_GLOB);
  ref_malloc(a_edge, a_total, REF_INT);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_total += b_size[part];
  ref_malloc(b_nodes, 2 * b_total, REF_GLOB);
  ref_malloc(b_data, b_total, REF_GLOB);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) a_next[part] =
      a_next[part - 1] + a_size[part - 1];

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (part != ref_mpi_rank(ref_mpi)) {
      a_edge[a_next[part]] = edge;
      a_nodes[0 + 2 * a_next[part]] =
          ref_node_global(ref_node, ref_edge_e2n(ref_edge, 0, edge));
      a_nodes[1 + 2 * a_next[part]] =
          ref_node_global(ref_node, ref_edge_e2n(ref_edge, 1, edge));
      (a_next[part])++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_nodes, a_size, b_nodes, b_size, 2,
                        REF_GLOB_TYPE),
      "alltoallv requested nodes");

  for (request = 0; request < b_total; request++) {
    RSS(ref_node_local(ref_node, b_nodes[0 + 2 * request], &node0), "loc 0");
    RSS(ref_node_local(ref_node, b_nodes[1 + 2 * request], &node1), "loc 1");
    RSS(ref_edge_with(ref_edge, node0, node1, &edge), "find edge");
    b_data[request] = data[edge];
  }

  RSS(ref_mpi_alltoallv(ref_mpi, b_data, b_size, a_data, a_size, 1,
                        REF_GLOB_TYPE),
      "alltoallv return data");

  for (request = 0; request < a_total; request++) {
    data[a_edge[request]] = a_data[request];
  }

  ref_free(a_next);

  ref_free(b_data);
  ref_free(b_nodes);

  ref_free(a_edge);

  ref_free(a_data);
  ref_free(a_nodes);

  ref_free(b_size);
  ref_free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_edge_ghost_min_int(REF_EDGE ref_edge, REF_MPI ref_mpi,
                                  REF_INT *data) {
  REF_NODE ref_node = ref_edge_node(ref_edge);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT edge;
  REF_INT part;
  REF_INT *a_next, *a_edge;
  REF_GLOB *a_nodes, *b_nodes;
  REF_INT *a_data, *b_data;

  REF_INT node0, node1;
  REF_INT request;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (part != ref_mpi_rank(ref_mpi)) a_size[part]++;
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_total += a_size[part];
  ref_malloc(a_nodes, 2 * a_total, REF_GLOB);
  ref_malloc(a_data, a_total, REF_INT);
  ref_malloc(a_edge, a_total, REF_INT);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_total += b_size[part];
  ref_malloc(b_nodes, 2 * b_total, REF_GLOB);
  ref_malloc(b_data, b_total, REF_INT);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) a_next[part] =
      a_next[part - 1] + a_size[part - 1];

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (part != ref_mpi_rank(ref_mpi)) {
      a_edge[a_next[part]] = edge;
      a_data[a_next[part]] = data[edge];
      a_nodes[0 + 2 * a_next[part]] =
          ref_node_global(ref_node, ref_edge_e2n(ref_edge, 0, edge));
      a_nodes[1 + 2 * a_next[part]] =
          ref_node_global(ref_node, ref_edge_e2n(ref_edge, 1, edge));
      (a_next[part])++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_nodes, a_size, b_nodes, b_size, 2,
                        REF_GLOB_TYPE),
      "alltoallv requested nodes");
  RSS(ref_mpi_alltoallv(ref_mpi, a_data, a_size, b_data, b_size, 1,
                        REF_INT_TYPE),
      "alltoallv requested data");

  /* min local data with remote ghost data */
  for (request = 0; request < b_total; request++) {
    RSS(ref_node_local(ref_node, b_nodes[0 + 2 * request], &node0), "loc 0");
    RSS(ref_node_local(ref_node, b_nodes[1 + 2 * request], &node1), "loc 1");
    RSS(ref_edge_with(ref_edge, node0, node1, &edge), "find edge");
    data[edge] = MIN(b_data[request], data[edge]);
  }
  /* export local consistent data */
  for (request = 0; request < b_total; request++) {
    RSS(ref_node_local(ref_node, b_nodes[0 + 2 * request], &node0), "loc 0");
    RSS(ref_node_local(ref_node, b_nodes[1 + 2 * request], &node1), "loc 1");
    RSS(ref_edge_with(ref_edge, node0, node1, &edge), "find edge");
    b_data[request] = data[edge];
  }

  RSS(ref_mpi_alltoallv(ref_mpi, b_data, b_size, a_data, a_size, 1,
                        REF_INT_TYPE),
      "alltoallv return data");

  for (request = 0; request < a_total; request++) {
    data[a_edge[request]] = a_data[request];
  }

  ref_free(a_next);

  ref_free(b_data);
  ref_free(b_nodes);

  ref_free(a_edge);

  ref_free(a_data);
  ref_free(a_nodes);

  ref_free(b_size);
  ref_free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_edge_ghost_dbl(REF_EDGE ref_edge, REF_MPI ref_mpi, REF_DBL *data,
                              REF_INT dim) {
  REF_NODE ref_node = ref_edge_node(ref_edge);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT edge;
  REF_INT part;
  REF_INT *a_next, *a_edge;
  REF_GLOB *a_nodes, *b_nodes;
  REF_DBL *a_data, *b_data;

  REF_INT node0, node1;
  REF_INT request;
  REF_INT i;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (part != ref_mpi_rank(ref_mpi)) a_size[part]++;
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_total += a_size[part];
  ref_malloc(a_nodes, 2 * a_total, REF_GLOB);
  ref_malloc(a_data, dim * a_total, REF_DBL);
  ref_malloc(a_edge, a_total, REF_INT);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_total += b_size[part];
  ref_malloc(b_nodes, 2 * b_total, REF_GLOB);
  ref_malloc(b_data, dim * b_total, REF_DBL);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) a_next[part] =
      a_next[part - 1] + a_size[part - 1];

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_edge_part(ref_edge, edge, &part), "edge part");
    if (part != ref_mpi_rank(ref_mpi)) {
      a_edge[a_next[part]] = edge;
      a_nodes[0 + 2 * a_next[part]] =
          ref_node_global(ref_node, ref_edge_e2n(ref_edge, 0, edge));
      a_nodes[1 + 2 * a_next[part]] =
          ref_node_global(ref_node, ref_edge_e2n(ref_edge, 1, edge));
      (a_next[part])++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_nodes, a_size, b_nodes, b_size, 2,
                        REF_GLOB_TYPE),
      "alltoallv requested nodes");

  for (request = 0; request < b_total; request++) {
    RSS(ref_node_local(ref_node, b_nodes[0 + 2 * request], &node0), "loc 0");
    RSS(ref_node_local(ref_node, b_nodes[1 + 2 * request], &node1), "loc 1");
    RSS(ref_edge_with(ref_edge, node0, node1, &edge), "find edge");
    for (i = 0; i < dim; i++) b_data[i + dim * request] = data[i + dim * edge];
  }

  RSS(ref_mpi_alltoallv(ref_mpi, b_data, b_size, a_data, a_size, dim,
                        REF_DBL_TYPE),
      "alltoallv return data");

  for (request = 0; request < a_total; request++) {
    for (i = 0; i < dim; i++)
      data[i + dim * a_edge[request]] = a_data[i + dim * request];
  }

  ref_free(a_next);

  ref_free(b_data);
  ref_free(b_nodes);

  ref_free(a_edge);

  ref_free(a_data);
  ref_free(a_nodes);

  ref_free(b_size);
  ref_free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_edge_tec_fill(REF_EDGE ref_edge, const char *filename) {
  REF_INT edge;

  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine edge fill\"\n");
  fprintf(file, "variables = \"i\" \"j\"\n");

  fprintf(file, "zone t=\"fill\", i=%d, datapacking=%s\n",
          2 * ref_edge_n(ref_edge), "point");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    fprintf(file, " %d %d\n", ref_edge_e2n(ref_edge, 0, edge),
            ref_edge_e2n(ref_edge, 1, edge));
    fprintf(file, " %d %d\n", ref_edge_e2n(ref_edge, 1, edge),
            ref_edge_e2n(ref_edge, 0, edge));
  }

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_edge_tec_int(REF_EDGE ref_edge, const char *filename,
                            REF_INT *data) {
  REF_NODE ref_node = ref_edge_node(ref_edge);
  REF_INT edge;
  REF_INT node;

  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine scalar file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\" \"s\"\n");

  fprintf(
      file,
      "zone t=\"scalar\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      2 * ref_edge_n(ref_edge), ref_edge_n(ref_edge), "point", "felineseg");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    node = ref_edge_e2n(ref_edge, 0, edge);
    fprintf(file, " %.16e %.16e %.16e %d\n", ref_node_xyz(ref_node, 0, node),
            ref_node_xyz(ref_node, 1, node), ref_node_xyz(ref_node, 2, node),
            data[edge]);
    node = ref_edge_e2n(ref_edge, 1, edge);
    fprintf(file, " %.16e %.16e %.16e %d\n", ref_node_xyz(ref_node, 0, node),
            ref_node_xyz(ref_node, 1, node), ref_node_xyz(ref_node, 2, node),
            data[edge]);
  }

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++)
    fprintf(file, " %d %d\n", 1 + 2 * edge, 2 + 2 * edge);

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_edge_tec_dbl(REF_EDGE ref_edge, const char *filename,
                            REF_DBL *data) {
  REF_NODE ref_node = ref_edge_node(ref_edge);
  REF_INT edge;
  REF_INT node;

  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine scalar file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\" \"r\"\n");

  fprintf(
      file,
      "zone t=\"scalar\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      2 * ref_edge_n(ref_edge), ref_edge_n(ref_edge), "point", "felineseg");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    node = ref_edge_e2n(ref_edge, 0, edge);
    fprintf(file, " %.16e %.16e %.16e %.16e\n", ref_node_xyz(ref_node, 0, node),
            ref_node_xyz(ref_node, 1, node), ref_node_xyz(ref_node, 2, node),
            data[edge]);
    node = ref_edge_e2n(ref_edge, 1, edge);
    fprintf(file, " %.16e %.16e %.16e %.16e\n", ref_node_xyz(ref_node, 0, node),
            ref_node_xyz(ref_node, 1, node), ref_node_xyz(ref_node, 2, node),
            data[edge]);
  }

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++)
    fprintf(file, " %d %d\n", 1 + 2 * edge, 2 + 2 * edge);

  fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_edge_tec_ratio(REF_EDGE ref_edge, const char *filename) {
  REF_NODE ref_node = ref_edge_node(ref_edge);
  REF_INT edge;
  REF_INT node;
  REF_DBL ratio;

  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine scalar file\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\" \"ratio\"\n");

  fprintf(
      file,
      "zone t=\"scalar\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      2 * ref_edge_n(ref_edge), ref_edge_n(ref_edge), "point", "felineseg");

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++) {
    RSS(ref_node_ratio(ref_node, ref_edge_e2n(ref_edge, 0, edge),
                       ref_edge_e2n(ref_edge, 1, edge), &ratio),
        "rat");
    node = ref_edge_e2n(ref_edge, 0, edge);
    fprintf(file, " %.16e %.16e %.16e %.16e\n", ref_node_xyz(ref_node, 0, node),
            ref_node_xyz(ref_node, 1, node), ref_node_xyz(ref_node, 2, node),
            ratio);
    node = ref_edge_e2n(ref_edge, 1, edge);
    fprintf(file, " %.16e %.16e %.16e %.16e\n", ref_node_xyz(ref_node, 0, node),
            ref_node_xyz(ref_node, 1, node), ref_node_xyz(ref_node, 2, node),
            ratio);
  }

  for (edge = 0; edge < ref_edge_n(ref_edge); edge++)
    fprintf(file, " %d %d\n", 1 + 2 * edge, 2 + 2 * edge);

  fclose(file);

  return REF_SUCCESS;
}

static REF_STATUS ref_edge_min_degree_node(REF_EDGE ref_edge, REF_INT *o2n,
                                           REF_INT *min_degree,
                                           REF_INT *min_degree_node) {
  REF_NODE ref_node = ref_edge_node(ref_edge);
  REF_ADJ ref_adj = ref_edge_adj(ref_edge);
  REF_INT node, degree;
  *min_degree = REF_EMPTY;
  *min_degree_node = REF_EMPTY;

  each_ref_node_valid_node(ref_node, node) {
    if (REF_EMPTY != o2n[node]) continue;
    RSS(ref_adj_degree(ref_adj, node, &degree), "deg");
    if (degree > 0) {
      if (REF_EMPTY == (*min_degree_node) || degree < (*min_degree)) {
        *min_degree_node = node;
        *min_degree = degree;
      }
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_edge_rcm_queue_node(REF_INT node, REF_INT degree,
                                          REF_INT *queue, REF_INT *nqueue,
                                          REF_INT *nhere) {
  REF_INT location, insert_point;

  /* largest degree first, will dequeue smallest from end of array */

  insert_point = *nhere;
  for (location = 0; location < (*nhere); location++) {
    if (queue[1 + 2 * location] < degree) {
      insert_point = location;
      break;
    }
  }

  for (location = (*nqueue); location > insert_point; location--) {
    queue[0 + 2 * location] = queue[0 + 2 * (location - 1)];
    queue[1 + 2 * location] = queue[1 + 2 * (location - 1)];
  }

  queue[0 + 2 * insert_point] = node;
  queue[1 + 2 * insert_point] = degree;
  (*nqueue)++;

  (*nhere)++;

  /*
  printf("nqueue %d of %d (%d) ins %d\n", *nqueue, node, degree, insert_point);
  for (location = 0; location < (*nqueue); location++) {
    printf("%d: %d (%d)\n", location, queue[0 + 2 * location],
           queue[1 + 2 * location]);
  }
  */

  return REF_SUCCESS;
}

static REF_STATUS ref_edge_rcm_queue(REF_EDGE ref_edge, REF_INT node,
                                     REF_INT *o2n, REF_INT *n2o, REF_INT *ndone,
                                     REF_INT *queue, REF_INT *nqueue) {
  REF_ADJ ref_adj = ref_edge_adj(ref_edge);
  REF_INT item, ref, other, degree;
  REF_INT nhere;
  n2o[(*ndone)] = node;
  o2n[node] = (*ndone);
  (*ndone)++;
  nhere = 0;

  /* printf("%d done %d nq %d\n", *ndone, node, *nqueue); */

  each_ref_adj_node_item_with_ref(ref_adj, node, item, ref) {
    other = ref_edge_e2n(ref_edge, 0, ref);
    if (REF_EMPTY == o2n[other]) {
      o2n[other] = -2; /* mark as queued */
      RSS(ref_adj_degree(ref_adj, other, &degree), "deg");
      RSS(ref_edge_rcm_queue_node(other, degree, queue, nqueue, &nhere),
          "queue n0");
    }
    other = ref_edge_e2n(ref_edge, 1, ref);
    if (REF_EMPTY == o2n[other]) {
      o2n[other] = -2; /* mark as queued */
      RSS(ref_adj_degree(ref_adj, other, &degree), "deg");
      RSS(ref_edge_rcm_queue_node(other, degree, queue, nqueue, &nhere),
          "queue n1");
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_edge_rcm(REF_EDGE ref_edge, REF_INT **o2n_ptr,
                        REF_INT **n2o_ptr) {
  REF_NODE ref_node = ref_edge_node(ref_edge);
  REF_INT *o2n, *n2o, *queue;
  REF_INT min_degree, min_degree_node;
  REF_INT ndone, nqueue, node;

  *o2n_ptr = NULL;
  *n2o_ptr = NULL;
  ref_malloc_init(o2n, ref_node_max(ref_node), REF_INT, REF_EMPTY);
  ref_malloc(n2o, ref_node_n(ref_node), REF_INT);
  ref_malloc(queue, 2 * ref_node_n(ref_node), REF_INT);

  ndone = 0;
  nqueue = 0;

  while (ndone < ref_node_n(ref_node)) {
    RSS(ref_edge_min_degree_node(ref_edge, o2n, &min_degree, &min_degree_node),
        "min degree node");
    RSS(ref_edge_rcm_queue(ref_edge, min_degree_node, o2n, n2o, &ndone, queue,
                           &nqueue),
        "min");

    /* drain queue */
    while (nqueue > 0) {
      min_degree_node = queue[0 + 2 * (nqueue - 1)];
      min_degree = queue[1 + 2 * (nqueue - 1)];
      nqueue--;
      RSS(ref_edge_rcm_queue(ref_edge, min_degree_node, o2n, n2o, &ndone, queue,
                             &nqueue),
          "min");
    }
  }

  REIS(ndone, ref_node_n(ref_node), "reordering done not original nodes");

  /* reverse with queue as temporary space */
  for (node = 0; node < ndone; node++) {
    queue[ndone - node - 1] = n2o[node];
  }
  /* copy back */
  for (node = 0; node < ndone; node++) {
    n2o[node] = queue[node];
  }
  ref_free(queue);

  /* set o2n to reverse */
  for (node = 0; node < ndone; node++) {
    o2n[n2o[node]] = node;
  }

  *o2n_ptr = o2n;
  *n2o_ptr = n2o;

  return REF_SUCCESS;
}
