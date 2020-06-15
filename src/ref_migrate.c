
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
#undef HAVE_MPI /* sometines defined by zoltan.h */
#include "zoltan.h"
#ifndef HAVE_MPI
#define HAVE_MPI
#endif
#endif

#if defined(HAVE_PARMETIS) && defined(HAVE_MPI)
#include "mpi.h"
#include "parmetis.h"
#if PARMETIS_MAJOR_VERSION == 3
#define PARM_INT idxtype
#define PARM_REAL float
#else
#define PARM_INT idx_t
#define PARM_REAL real_t
#endif
#endif

#include "ref_export.h"
#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_migrate.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_sort.h"

REF_STATUS ref_migrate_create(REF_MIGRATE *ref_migrate_ptr, REF_GRID ref_grid) {
  REF_MIGRATE ref_migrate;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT node;
  REF_INT group, cell, cell_edge, n0, n1;

  ref_malloc(*ref_migrate_ptr, 1, REF_MIGRATE_STRUCT);

  ref_migrate = *ref_migrate_ptr;

  ref_migrate_grid(ref_migrate) = ref_grid;

  RSS(ref_adj_create(&(ref_migrate_parent_local(ref_migrate))), "make adj");
  RSS(ref_adj_create(&(ref_migrate_parent_part(ref_migrate))), "make adj");
  RSS(ref_adj_create(&(ref_migrate_conn(ref_migrate))), "make adj");

  ref_migrate_max(ref_migrate) = ref_node_max(ref_node);

  ref_malloc_init(ref_migrate->global, ref_migrate_max(ref_migrate), REF_GLOB,
                  REF_EMPTY);
  ref_malloc(ref_migrate->xyz, 3 * ref_migrate_max(ref_migrate), REF_DBL);
  ref_malloc(ref_migrate->weight, ref_migrate_max(ref_migrate), REF_DBL);
  ref_malloc(ref_migrate->age, ref_migrate_max(ref_migrate), REF_INT);

  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      ref_migrate_global(ref_migrate, node) = ref_node_global(ref_node, node);
      RSS(ref_adj_add(ref_migrate_parent_local(ref_migrate), node, node),
          "add");
      RSS(ref_adj_add(ref_migrate_parent_part(ref_migrate), node,
                      ref_node_part(ref_node, node)),
          "add");
      ref_migrate_xyz(ref_migrate, 0, node) = ref_node_xyz(ref_node, 0, node);
      ref_migrate_xyz(ref_migrate, 1, node) = ref_node_xyz(ref_node, 1, node);
      ref_migrate_xyz(ref_migrate, 2, node) = ref_node_xyz(ref_node, 2, node);
      ref_migrate_weight(ref_migrate, node) = 1.0;
      ref_migrate_age(ref_migrate, node) = ref_node_age(ref_node, node);
    }
  }
  RSS(ref_node_ghost_int(ref_node, (ref_migrate->age), 1),
      "ghost age for edge weights");

  /* 2d included for twod */
  each_ref_grid_2d_3d_ref_cell(ref_grid, group, ref_cell) {
    each_ref_cell_valid_cell(ref_cell, cell) {
      each_ref_cell_cell_edge(ref_cell, cell_edge) {
        /* need ghost nodes for agglomeration */
        n0 = ref_cell_e2n(ref_cell, 0, cell_edge, cell);
        n1 = ref_cell_e2n(ref_cell, 1, cell_edge, cell);
        RSS(ref_adj_add_uniquely(ref_migrate_conn(ref_migrate), n0, n1),
            "uniq");
        RSS(ref_adj_add_uniquely(ref_migrate_conn(ref_migrate), n1, n0),
            "uniq");
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_free(REF_MIGRATE ref_migrate) {
  if (NULL == (void *)ref_migrate) return REF_NULL;

  ref_free(ref_migrate->age);
  ref_free(ref_migrate->weight);
  ref_free(ref_migrate->xyz);
  ref_free(ref_migrate->global);

  RSS(ref_adj_free(ref_migrate_conn(ref_migrate)), "free adj");
  RSS(ref_adj_free(ref_migrate_parent_part(ref_migrate)), "free adj");
  RSS(ref_adj_free(ref_migrate_parent_local(ref_migrate)), "free adj");

  ref_free(ref_migrate);

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_inspect(REF_MIGRATE ref_migrate) {
  REF_NODE ref_node = ref_grid_node(ref_migrate_grid(ref_migrate));
  REF_INT node, item, local, part;
  REF_GLOB global;

  each_ref_migrate_node(ref_migrate, node) {
    printf(" %2d : " REF_GLOB_FMT " :", ref_mpi_rank(ref_node_mpi(ref_node)),
           ref_node_global(ref_node, node));
    each_ref_adj_node_item_with_ref(ref_migrate_parent_local(ref_migrate), node,
                                    item, local) {
      global = ref_migrate_global(ref_migrate, local);
      part = ref_adj_item_ref(ref_migrate_parent_part(ref_migrate), item);
      printf(" " REF_GLOB_FMT "+%d", global, part);
    }
    printf("\n");
  }
  return REF_SUCCESS;
}

REF_STATUS ref_migrate_2d_agglomeration_keep(REF_MIGRATE ref_migrate,
                                             REF_INT keep, REF_INT lose) {
  REF_NODE ref_node = ref_grid_node(ref_migrate_grid(ref_migrate));
  REF_ADJ conn_adj = ref_migrate_conn(ref_migrate);
  REF_INT item, local;
  REF_GLOB global;
  REF_INT from_node;

  /* not working for general agglomeration, ghost lose? */

  RAS(ref_node_valid(ref_node, keep), "keep node invalid");
  RAS(ref_node_valid(ref_node, lose), "lose node invalid");

  ref_migrate_global(ref_migrate, lose) = REF_EMPTY;

  /* skip if the lose node has been agglomerated */
  each_ref_adj_node_item_with_ref(ref_migrate_parent_local(ref_migrate), keep,
                                  item, local) {
    global = ref_migrate_global(ref_migrate, local);
    if (global == ref_node_global(ref_node, lose)) {
      return REF_SUCCESS;
    }
  }
  /* update edges pointing to lose node */
  each_ref_adj_node_item_with_ref(conn_adj, lose, item, from_node) {
    RSS(ref_adj_remove(conn_adj, from_node, lose), "rm to lose");
    if (from_node != keep) {
      RSS(ref_adj_add_uniquely(conn_adj, from_node, keep), "add to keep");
      RSS(ref_adj_add_uniquely(conn_adj, keep, from_node), "add to keep");
    }
  }

  /* update edges pointing from lose node */
  while (ref_adj_valid(ref_adj_first(conn_adj, lose))) {
    RSS(ref_adj_remove(
            conn_adj, lose,
            ref_adj_item_ref(conn_adj, ref_adj_first(conn_adj, lose))),
        "rm from lose");
  }

  /* skip if keep node is off-proc or already agglomerated */
  if (!ref_migrate_valid(ref_migrate, keep)) return REF_SUCCESS;

  ref_migrate_xyz(ref_migrate, 1, keep) = 0.5;
  ref_migrate_weight(ref_migrate, keep) = 2.0;
  /* collect age in general case */
  RSS(ref_adj_add(ref_migrate_parent_local(ref_migrate), keep, lose), "add");
  RSS(ref_adj_add(ref_migrate_parent_part(ref_migrate), keep,
                  ref_node_part(ref_node, lose)),
      "add");

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_2d_agglomeration(REF_MIGRATE ref_migrate) {
  REF_GRID ref_grid = ref_migrate_grid(ref_migrate);
  REF_NODE ref_node = ref_grid_node(ref_migrate_grid(ref_migrate));
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT keep, lose;

  each_ref_cell_valid_cell_with_nodes(ref_grid_pri(ref_grid), cell, nodes) {
    if (ref_node_global(ref_node, nodes[0]) <
        ref_node_global(ref_node, nodes[3])) {
      keep = nodes[0];
      lose = nodes[3];
    } else {
      keep = nodes[3];
      lose = nodes[0];
    }
    RSS(ref_migrate_2d_agglomeration_keep(ref_migrate, keep, lose), "0-3");

    if (ref_node_global(ref_node, nodes[1]) <
        ref_node_global(ref_node, nodes[4])) {
      keep = nodes[1];
      lose = nodes[4];
    } else {
      keep = nodes[4];
      lose = nodes[1];
    }
    RSS(ref_migrate_2d_agglomeration_keep(ref_migrate, keep, lose), "1-4");

    if (ref_node_global(ref_node, nodes[2]) <
        ref_node_global(ref_node, nodes[5])) {
      keep = nodes[2];
      lose = nodes[5];
    } else {
      keep = nodes[5];
      lose = nodes[2];
    }
    RSS(ref_migrate_2d_agglomeration_keep(ref_migrate, keep, lose), "2-5");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_migrate_report_load_balance(REF_GRID ref_grid,
                                                  REF_INT *node_part) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_INT min_part, max_part, node, proc, *partition_size;
  ref_malloc_init(partition_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      RAB(0 <= node_part[node] && node_part[node] < ref_mpi_n(ref_mpi),
          "part out of range", {
            printf("rank %d node %d node_part %d n %d", ref_mpi_rank(ref_mpi),
                   node, node_part[node], ref_mpi_n(ref_mpi));
          });
      partition_size[node_part[node]] += 1;
    }
  }
  RSS(ref_mpi_allsum(ref_mpi, partition_size, ref_mpi_n(ref_mpi), REF_INT_TYPE),
      "allsum");

  min_part = INT_MAX;
  max_part = 0;
  each_ref_mpi_part(ref_mpi, proc) {
    min_part = MIN(min_part, partition_size[proc]);
    max_part = MAX(max_part, partition_size[proc]);
  }

  if (ref_mpi_once(ref_mpi)) {
    printf(
        "balance %6.3f on %d target %d size min %d max %d\n",
        (REF_DBL)max_part / (REF_DBL)ref_node_n_global(ref_node) *
            (REF_DBL)ref_mpi_n(ref_mpi),
        ref_mpi_n(ref_mpi),
        (REF_INT)(ref_node_n_global(ref_node) / (REF_GLOB)ref_mpi_n(ref_mpi)),
        min_part, max_part);
  }
  ref_free(partition_size);
  return REF_SUCCESS;
}

static REF_STATUS ref_migrate_single_part(REF_GRID ref_grid,
                                          REF_INT *node_part) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;

  for (node = 0; node < ref_node_max(ref_node); node++) node_part[node] = 0;

  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "single part");
  return REF_SUCCESS;
}

static REF_STATUS ref_migrate_native_rcb_direction(
    REF_MPI ref_mpi, REF_INT n, REF_DBL *xyz, REF_INT npart, REF_INT *owners,
    REF_INT *locals, REF_MPI global_mpi, REF_INT *part, REF_INT seed,
    REF_INT dir) {
  REF_INT i, j, n0, n1, npart0, npart1;
  REF_INT bal_n0, bal_n1;
  REF_DBL *xyz0, *xyz1, *x;
  REF_DBL *bal_xyz0, *bal_xyz1;
  REF_INT *owners0, *owners1;
  REF_INT *bal_owners0, *bal_owners1;
  REF_INT *locals0, *locals1;
  REF_INT *bal_locals0, *bal_locals1;
  REF_DBL ratio, value0, value1;
  REF_LONG position, total;
  REF_MPI split_mpi;
  REF_INT seed_base = 3;
  REF_DBL ratio_shift, ratio0, ratio1;

  if (0 == npart) return REF_SUCCESS;

  if (1 == npart) {
    REF_INT *my_id, *recv_part, *recv_locals, nrecv;
    ref_malloc_init(my_id, n, REF_INT, ref_mpi_rank(global_mpi));
    RSS(ref_mpi_blindsend(global_mpi, owners, my_id, 1, n, (void **)&recv_part,
                          &nrecv, REF_INT_TYPE),
        "recv part");
    RSS(ref_mpi_blindsend(global_mpi, owners, locals, 1, n,
                          (void **)&recv_locals, &nrecv, REF_INT_TYPE),
        "recv loc");
    for (i = 0; i < nrecv; i++) part[recv_locals[i]] = recv_part[i];
    ref_free(recv_part);
    ref_free(recv_locals);
    ref_free(my_id);
    return REF_SUCCESS;
  }

  ref_malloc(x, n, REF_DBL);
  if (dir < 0 || 2 < dir)
    RSS(ref_migrate_split_dir(ref_mpi, n, xyz, &dir), "dir");
  RAS(-1 < dir && dir < 3, "3D dir");
  RSS(ref_migrate_split_ratio(npart, &ratio), "ratio");
  ratio_shift = (REF_DBL)(seed % seed_base) / (REF_DBL)seed_base;
  ratio0 = ratio * ratio_shift;
  ratio1 = 1.0 - (ratio - ratio0);

  for (i = 0; i < n; i++) x[i] = xyz[dir + 3 * i];

  total = (REF_LONG)n;
  RSS(ref_mpi_allsum(ref_mpi, &total, 1, REF_LONG_TYPE), "high_pos");

  position = (REF_LONG)((REF_DBL)total * ratio0);
  RSS(ref_search_selection(ref_mpi, n, x, position, &value0), "target");
  position = (REF_LONG)((REF_DBL)total * ratio1);
  RSS(ref_search_selection(ref_mpi, n, x, position, &value1), "target");

  ref_malloc(xyz0, 3 * n, REF_DBL);
  ref_malloc(xyz1, 3 * n, REF_DBL);
  ref_malloc(owners0, n, REF_INT);
  ref_malloc(owners1, n, REF_INT);
  ref_malloc(locals0, n, REF_INT);
  ref_malloc(locals1, n, REF_INT);

  n0 = 0;
  n1 = 0;
  for (i = 0; i < n; i++) {
    if (x[i] < value0 || value1 < x[i]) {
      for (j = 0; j < 3; j++) xyz0[j + 3 * n0] = xyz[j + 3 * i];
      owners0[n0] = owners[i];
      locals0[n0] = locals[i];
      n0++;
    } else {
      for (j = 0; j < 3; j++) xyz1[j + 3 * n1] = xyz[j + 3 * i];
      owners1[n1] = owners[i];
      locals1[n1] = locals[i];
      n1++;
    }
  }
  REIS(n, n0 + n1, "conservation");
  npart0 = npart / 2;
  npart1 = npart - npart0;

  RSS(ref_mpi_balance(ref_mpi, 3, n0, (void *)xyz0, 0, npart0 - 1, &bal_n0,
                      (void **)(&bal_xyz0), REF_DBL_TYPE),
      "split 0");
  RSS(ref_mpi_balance(ref_mpi, 3, n1, (void *)xyz1, npart0,
                      ref_mpi_n(ref_mpi) - 1, &bal_n1, (void **)(&bal_xyz1),
                      REF_DBL_TYPE),
      "split 1");

  RSS(ref_mpi_balance(ref_mpi, 1, n0, (void *)owners0, 0, npart0 - 1, &bal_n0,
                      (void **)(&bal_owners0), REF_INT_TYPE),
      "split owner 0");
  RSS(ref_mpi_balance(ref_mpi, 1, n1, (void *)owners1, npart0,
                      ref_mpi_n(ref_mpi) - 1, &bal_n1, (void **)(&bal_owners1),
                      REF_INT_TYPE),
      "split owner 1");

  RSS(ref_mpi_balance(ref_mpi, 1, n0, (void *)locals0, 0, npart0 - 1, &bal_n0,
                      (void **)(&bal_locals0), REF_INT_TYPE),
      "split local 0");
  RSS(ref_mpi_balance(ref_mpi, 1, n1, (void *)locals1, npart0,
                      ref_mpi_n(ref_mpi) - 1, &bal_n1, (void **)(&bal_locals1),
                      REF_INT_TYPE),
      "split local 1");

  RSS(ref_mpi_front_comm(ref_mpi, &split_mpi, npart0), "split");

  dir += 1;
  if (dir > 2) dir -= 3;
  if (ref_mpi_rank(ref_mpi) < npart0) {
    RSS(ref_migrate_native_rcb_direction(split_mpi, bal_n0, bal_xyz0, npart0,
                                         bal_owners0, bal_locals0, global_mpi,
                                         part, seed, dir),
        "recurse 0");
  } else {
    RSS(ref_migrate_native_rcb_direction(split_mpi, bal_n1, bal_xyz1, npart1,
                                         bal_owners1, bal_locals1, global_mpi,
                                         part, seed, dir),
        "recurse 1");
  }

  RSS(ref_mpi_join_comm(split_mpi), "join");
  RSS(ref_mpi_free(split_mpi), "new free");

  ref_free(bal_locals1);
  ref_free(bal_locals0);

  ref_free(bal_owners1);
  ref_free(bal_owners0);

  ref_free(bal_xyz1);
  ref_free(bal_xyz0);

  ref_free(locals1);
  ref_free(locals0);

  ref_free(owners1);
  ref_free(owners0);

  ref_free(xyz1);
  ref_free(xyz0);

  ref_free(x);
  return REF_SUCCESS;
}

static REF_STATUS ref_migrate_native_rcb_part(REF_GRID ref_grid,
                                              REF_INT *node_part) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_INT node;
  REF_INT i, n;
  REF_DBL *xyz;
  REF_INT npart;
  REF_INT *owners;
  REF_INT *locals;

  for (node = 0; node < ref_node_max(ref_node); node++)
    node_part[node] = REF_EMPTY;

  npart = ref_mpi_n(ref_mpi);

  n = ref_node_n(ref_node);
  ref_malloc(xyz, 3 * n, REF_DBL);
  ref_malloc(owners, n, REF_INT);
  ref_malloc(locals, n, REF_INT);
  n = 0;
  each_ref_node_valid_node(ref_node, node) {
    if (ref_node_owned(ref_node, node)) {
      for (i = 0; i < 3; i++) xyz[i + 3 * n] = ref_node_xyz(ref_node, i, node);
      owners[n] = ref_node_part(ref_node, node);
      locals[n] = node;
      n++;
    }
  }

  RSS(ref_migrate_native_rcb_direction(ref_mpi, n, xyz, npart, owners, locals,
                                       ref_mpi, node_part,
                                       ref_grid_partitioner_seed(ref_grid), -1),
      "split");
  ref_grid_partitioner_seed(ref_grid)++;
  if (ref_grid_partitioner_seed(ref_grid) < 0)
    ref_grid_partitioner_seed(ref_grid) = 0; /* overflow int */

  ref_free(locals);
  ref_free(owners);
  ref_free(xyz);

  RSS(ref_migrate_report_load_balance(ref_grid, node_part), "report bal");

  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "native RCB part");

  return REF_SUCCESS;
}

#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
static int ref_migrate_zoltan_local_n(void *void_ref_migrate, int *ierr) {
  REF_MIGRATE ref_migrate = ((REF_MIGRATE)void_ref_migrate);
  int node, n;

  *ierr = 0;

  n = 0;
  each_ref_migrate_node(ref_migrate, node) { n++; }

  return n;
}
static void ref_migrate_zoltan_local_ids(void *void_ref_migrate, int global_dim,
                                         int local_dim, ZOLTAN_ID_PTR global,
                                         ZOLTAN_ID_PTR local, int wgt_dim,
                                         float *obj_wgts, int *ierr) {
  REF_MIGRATE ref_migrate = ((REF_MIGRATE)void_ref_migrate);
  REF_INT node, n;

  if (1 != global_dim || 1 != local_dim || 1 != wgt_dim) {
    printf("%s: %d: %s: %s\n", __FILE__, __LINE__, __func__, "bad sizes");
    *ierr = ZOLTAN_FATAL;
    return;
  }

  *ierr = 0;

  n = 0;
  each_ref_migrate_node(ref_migrate, node) {
    local[n] = (ZOLTAN_ID_TYPE)node;
    global[n] = (ZOLTAN_ID_TYPE)ref_migrate_global(ref_migrate, node);
    obj_wgts[n] = (float)ref_migrate_weight(ref_migrate, node);
    n++;
  }
}
static int ref_migrate_zoltan_geom_dimensionality(void *void_ref_migrate,
                                                  int *ierr) {
  SUPRESS_UNUSED_COMPILER_WARNING(void_ref_migrate);
  *ierr = 0;
  return 3;
}
static void ref_migrate_zoltan_geom(void *void_ref_migrate, int global_dim,
                                    int local_dim, int nnode,
                                    ZOLTAN_ID_PTR global, ZOLTAN_ID_PTR local,
                                    int xyz_dim, double *xyz, int *ierr) {
  REF_MIGRATE ref_migrate = ((REF_MIGRATE)void_ref_migrate);
  REF_INT node;

  SUPRESS_UNUSED_COMPILER_WARNING(global);
  *ierr = 0;

  if (1 != global_dim || 1 != local_dim || 3 != xyz_dim) {
    printf("%s: %d: %s: %s\n", __FILE__, __LINE__, __func__, "bad sizes");
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (node = 0; node < nnode; node++) {
    if (!ref_migrate_valid(ref_migrate, local[node])) {
      printf("%s: %d: %s: %d %d invalid\n", __FILE__, __LINE__, __func__, node,
             (REF_INT)local[node]);
      *ierr = ZOLTAN_FATAL;
      return;
    }
    xyz[0 + 3 * node] = ref_migrate_xyz(ref_migrate, 0, local[node]);
    xyz[1 + 3 * node] = ref_migrate_xyz(ref_migrate, 1, local[node]);
    xyz[2 + 3 * node] = ref_migrate_xyz(ref_migrate, 2, local[node]);
  }
}
static int ref_migrate_zoltan_num_edges(void *void_ref_migrate, int global_dim,
                                        int local_dim, ZOLTAN_ID_PTR global,
                                        ZOLTAN_ID_PTR local, int *ierr) {
  REF_MIGRATE ref_migrate = ((REF_MIGRATE)void_ref_migrate);
  REF_INT node, degree;

  SUPRESS_UNUSED_COMPILER_WARNING(global);
  *ierr = 0;

  if (1 != global_dim || 1 != local_dim) {
    printf("%s: %d: %s: %s\n", __FILE__, __LINE__, __func__, "bad sizes");
    *ierr = ZOLTAN_FATAL;
    return 0;
  }

  node = (REF_INT)local[0];
  RSS(ref_adj_degree(ref_migrate_conn(ref_migrate), node, &degree), "deg");

  return degree;
}
static void ref_migrate_zoltan_edge_list(void *void_ref_migrate, int global_dim,
                                         int local_dim, ZOLTAN_ID_PTR global,
                                         ZOLTAN_ID_PTR local,
                                         ZOLTAN_ID_PTR conn_global,
                                         int *conn_part, int weight_dim,
                                         float *weight, int *ierr) {
  REF_MIGRATE ref_migrate = ((REF_MIGRATE)void_ref_migrate);
  REF_NODE ref_node = ref_grid_node(ref_migrate_grid(ref_migrate));
  REF_INT node, item, ref, degree;

  SUPRESS_UNUSED_COMPILER_WARNING(global);
  SUPRESS_UNUSED_COMPILER_WARNING(weight);
  *ierr = 0;

  if (1 != global_dim || 1 != local_dim || 1 != weight_dim) {
    printf("%s: %d: %s: %s\n", __FILE__, __LINE__, __func__, "bad sizes");
    *ierr = ZOLTAN_FATAL;
    return;
  }

  node = (REF_INT)local[0];
  degree = 0;

  each_ref_adj_node_item_with_ref(ref_migrate_conn(ref_migrate), node, item,
                                  ref) {
    conn_global[degree] = (ZOLTAN_ID_TYPE)ref_node_global(ref_node, ref);
    conn_part[degree] = (int)ref_node_part(ref_node, ref);
    weight[degree] = (float)ref_migrate_age(ref_migrate, node) +
                     (float)ref_migrate_age(ref_migrate, ref) + (float)1.0;
    degree++;
  }
}
REF_STATUS ref_migrate_zoltan_part(REF_GRID ref_grid, REF_INT *node_part) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_MIGRATE ref_migrate;
  int partitions_have_changed;
  int global_id_dimension, local_id_dimension;

  int import_n;
  ZOLTAN_ID_PTR import_global, import_local;
  int *import_proc, *import_part;

  int export_n;
  ZOLTAN_ID_PTR export_global, export_local;
  int *export_proc, *export_part;

  float ver;

  REF_INT node, item, local, part;
  REF_GLOB global;

  REF_INT *migrate_part;

  REF_INT *a_next;
  REF_GLOB *a_parts, *b_parts;
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;

  struct Zoltan_Struct *zz;

  RSS(ref_node_synchronize_globals(ref_node), "sync global nodes");
  RSS(ref_node_collect_ghost_age(ref_node), "collect ghost age");

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  RSS(ref_migrate_create(&ref_migrate, ref_grid), "create migrate");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "zoltan init");

  if (ref_grid_twod(ref_grid)) {
    RSS(ref_migrate_2d_agglomeration(ref_migrate), "2d agglom");
  }

  { /* zoltan does not use argc and argv when MPI_Initialized
       must be protected */
    char **empty_argument = NULL;
    if (!ref_mpi_para(ref_mpi))
      THROW("Zoltan_Initialize must have actual arguments for seq");
    REIS(ZOLTAN_OK, Zoltan_Initialize(0, empty_argument, &ver),
         "Zoltan is angry");
  }
  zz = Zoltan_Create((*((MPI_Comm *)(ref_mpi->comm))));

  /* General parameters */

  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTS");
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");

  switch (ref_grid_partitioner(ref_grid)) {
    case REF_MIGRATE_RECOMMENDED:
    case REF_MIGRATE_ZOLTAN_GRAPH:
      Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
      break;
    case REF_MIGRATE_ZOLTAN_RCB:
      Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
      break;
    default:
      RSS(REF_IMPLEMENT, "ref_migrate_method");
      break;
  }

  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");

  Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1");

  Zoltan_Set_Num_Obj_Fn(zz, ref_migrate_zoltan_local_n, (void *)ref_migrate);
  Zoltan_Set_Obj_List_Fn(zz, ref_migrate_zoltan_local_ids, (void *)ref_migrate);
  Zoltan_Set_Num_Geom_Fn(zz, ref_migrate_zoltan_geom_dimensionality,
                         (void *)ref_migrate);
  Zoltan_Set_Geom_Multi_Fn(zz, ref_migrate_zoltan_geom, (void *)ref_migrate);

  Zoltan_Set_Num_Edges_Fn(zz, ref_migrate_zoltan_num_edges,
                          (void *)ref_migrate);
  Zoltan_Set_Edge_List_Fn(zz, ref_migrate_zoltan_edge_list,
                          (void *)ref_migrate);

  REIS(ZOLTAN_OK,
       Zoltan_LB_Partition(zz, &partitions_have_changed, &global_id_dimension,
                           &local_id_dimension, &import_n, &import_global,
                           &import_local, &import_proc, &import_part, &export_n,
                           &export_global, &export_local, &export_proc,
                           &export_part),
       "Zoltan is angry");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "zoltan part");

  ref_malloc_init(migrate_part, ref_migrate_max(ref_node), REF_INT, REF_EMPTY);

  for (node = 0; node < export_n; node++)
    migrate_part[export_local[node]] = export_part[node];

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_migrate_node(ref_migrate, node) {
    each_ref_adj_node_item_with_ref(ref_migrate_parent_local(ref_migrate), node,
                                    item, local) {
      part = ref_adj_item_ref(ref_migrate_parent_part(ref_migrate), item);
      if (ref_mpi_rank(ref_mpi) != part) {
        a_size[part]++;
      } else {
        node_part[local] = migrate_part[node];
      }
    }
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_total += a_size[part];
  ref_malloc(a_parts, 2 * a_total, REF_GLOB);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_total += b_size[part];
  ref_malloc(b_parts, 2 * b_total, REF_GLOB);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) {
    a_next[part] = a_next[part - 1] + a_size[part - 1];
  }

  each_ref_migrate_node(ref_migrate, node) {
    each_ref_adj_node_item_with_ref(ref_migrate_parent_local(ref_migrate), node,
                                    item, local) {
      part = ref_adj_item_ref(ref_migrate_parent_part(ref_migrate), item);
      if (ref_mpi_rank(ref_mpi) != part) {
        global = ref_migrate_global(ref_migrate, local);
        a_parts[0 + 2 * a_next[part]] = global;
        a_parts[1 + 2 * a_next[part]] = (REF_GLOB)migrate_part[node];
        a_next[part]++;
      }
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_parts, a_size, b_parts, b_size, 2,
                        REF_GLOB_TYPE),
      "alltoallv parts");

  for (node = 0; node < b_total; node++) {
    global = b_parts[0 + 2 * node];
    part = (REF_INT)b_parts[1 + 2 * node];
    RSS(ref_node_local(ref_node, global, &local), "g2l");
    node_part[local] = part;
  }

  free(a_next);
  free(b_parts);
  free(a_parts);
  free(b_size);
  free(a_size);

  ref_free(migrate_part);

  REIS(ZOLTAN_OK,
       Zoltan_LB_Free_Part(&import_local, &import_global, &import_proc,
                           &import_part),
       "Zoltan is angry");

  REIS(ZOLTAN_OK,
       Zoltan_LB_Free_Part(&export_local, &export_global, &export_proc,
                           &export_part),
       "Zoltan is angry");

  Zoltan_Destroy(&zz);

  RSS(ref_migrate_free(ref_migrate), "free migrate");

  RSS(ref_migrate_report_load_balance(ref_grid, node_part), "report bal");

  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "part update");

  return REF_SUCCESS;
}
#endif

#if defined(HAVE_PARMETIS) && defined(HAVE_MPI)
static REF_STATUS ref_migrate_metis_wrapper(REF_MPI ref_mpi, PARM_INT *vtxdist,
                                            PARM_INT *xadjdist,
                                            PARM_INT *adjncydist,
                                            PARM_INT *adjwgtdist,
                                            PARM_INT *partdist) {
  REF_INT *count;
  PARM_INT n, *xadj, *adjncy, *adjwgt, *part;
  PARM_INT *vwgt, *vsize, nparts, ncon, objval;
  PARM_REAL *tpwgts, *ubvec;
  PARM_INT options[METIS_NOPTIONS];
  REF_INT i, proc;

  n = vtxdist[ref_mpi_n(ref_mpi)];
  ref_malloc_init(count, ref_mpi_n(ref_mpi), REF_INT, REF_EMPTY);
  ref_malloc_init(xadj, n + 1, REF_INT, REF_EMPTY);
  each_ref_mpi_part(ref_mpi, proc) {
    count[proc] = vtxdist[proc + 1] - vtxdist[proc];
  }
  RSS(ref_mpi_allgatherv(ref_mpi, &(xadjdist[1]), count, &(xadj[1]),
                         REF_INT_TYPE),
      "gather adj");
  xadj[0] = 0;
  each_ref_mpi_part(ref_mpi, proc) {
    for (i = vtxdist[proc] + 1; i <= vtxdist[proc + 1]; i++) {
      xadj[i] += xadj[vtxdist[proc]];
    }
  }
  ref_malloc_init(adjncy, xadj[n], REF_INT, REF_EMPTY);
  ref_malloc_init(adjwgt, xadj[n], REF_INT, REF_EMPTY);
  each_ref_mpi_part(ref_mpi, proc) {
    count[proc] = xadj[vtxdist[proc + 1]] - xadj[vtxdist[proc]];
  }
  RSS(ref_mpi_allgatherv(ref_mpi, adjncydist, count, adjncy, REF_INT_TYPE),
      "gather adjncy");
  RSS(ref_mpi_allgatherv(ref_mpi, adjwgtdist, count, adjwgt, REF_INT_TYPE),
      "gather adjwgt");

  ref_mpi_stopwatch_stop(ref_mpi, "metis gather");

  ref_malloc_init(part, n, PARM_INT, REF_EMPTY);

  ncon = 1;
  vsize = NULL;
  nparts = ref_mpi_n(ref_mpi);

  ref_malloc_init(vwgt, ncon * n, PARM_INT, 1);
  ref_malloc_init(tpwgts, ncon * ref_mpi_n(ref_mpi), PARM_REAL,
                  1.0 / (PARM_REAL)ref_mpi_n(ref_mpi));
  ref_malloc_init(ubvec, ncon, PARM_REAL, 1.001);

  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_SEED] = 42;
  options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB; /* zero part less likely */
  /* options[METIS_OPTION_DBGLVL] = METIS_DBG_COARSEN; */
  if (ref_mpi_once(ref_mpi)) {
    REIS(METIS_OK,
         METIS_PartGraphKway(&n, &ncon, xadj, adjncy, vwgt, vsize, adjwgt,
                             &nparts, tpwgts, ubvec, options, &objval, part),
         "METIS is not o.k.");
  }

  ref_free(ubvec);
  ref_free(tpwgts);
  ref_free(vwgt);

  each_ref_mpi_part(ref_mpi, proc) {
    count[proc] = (REF_INT)(vtxdist[proc + 1] - vtxdist[proc]);
  }
  if (ref_mpi_once(ref_mpi)) {
    for (i = 0; i < vtxdist[1]; i++) {
      partdist[i] = part[i];
    }
    each_ref_mpi_worker(ref_mpi, proc) {
      RSS(ref_mpi_send(ref_mpi, &(part[vtxdist[proc]]), count[proc],
                       REF_INT_TYPE, proc),
          "send part");
    }
  } else {
    proc = ref_mpi_rank(ref_mpi);
    RSS(ref_mpi_recv(ref_mpi, partdist, count[proc], REF_INT_TYPE, 0),
        "recv part");
  }

  ref_free(part);
  ref_free(adjwgt);
  ref_free(adjncy);
  ref_free(xadj);
  ref_free(count);

  return REF_SUCCESS;
}
static REF_STATUS ref_migrate_parmetis_wrapper(
    REF_MPI ref_mpi, PARM_INT *vtxdist, PARM_INT *xadjdist,
    PARM_INT *adjncydist, PARM_INT *adjwgtdist, PARM_INT *partdist) {
  PARM_INT *vwgt;
  PARM_REAL *tpwgts, *ubvec;
  PARM_INT wgtflag = 3;
  PARM_INT numflag = 0;
  PARM_INT ncon;
  PARM_INT nparts;
  PARM_INT edgecut;
  PARM_INT options[] = {1, 0 /* PARMETIS_DBGLVL_PROGRESS */, 42};
  MPI_Comm comm = (*((MPI_Comm *)(ref_mpi->comm)));
  REF_INT n, proc;

  nparts = ref_mpi_n(ref_mpi);
  proc = ref_mpi_rank(ref_mpi);
  n = (REF_INT)(vtxdist[proc + 1] - vtxdist[proc]);
  ncon = 1;
  ref_malloc_init(vwgt, ncon * n, PARM_INT, 1);
  ref_malloc_init(tpwgts, ncon * ref_mpi_n(ref_mpi), PARM_REAL,
                  1.0 / (PARM_REAL)ref_mpi_n(ref_mpi));
  ref_malloc_init(ubvec, ncon, PARM_REAL, 1.01);

  REIS(METIS_OK,
       ParMETIS_V3_PartKway(vtxdist, xadjdist, adjncydist, vwgt, adjwgtdist,
                            &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec,
                            options, &edgecut, partdist, &comm),
       "ParMETIS is not o.k.");

  ref_free(ubvec);
  ref_free(tpwgts);
  ref_free(vwgt);
  return REF_SUCCESS;
}
static REF_STATUS ref_migrate_parmetis_subset(
    REF_MPI ref_mpi, REF_INT newproc, PARM_INT *vtxdist, PARM_INT *xadjdist,
    PARM_INT *adjncydist, PARM_INT *adjwgtdist, PARM_INT *partdist) {
  REF_INT proc, nold, nnew, i, first;
  REF_INT nsend, nrecv, *send_size, *recv_size;
  PARM_INT ntotal;
  PARM_INT n0, n1;
  PARM_INT *vtx, *xadj, *adjncy, *adjwgt, *part;
  PARM_INT *deg, *newdeg;
  REF_MPI split_mpi;
  ntotal = vtxdist[ref_mpi_n(ref_mpi)];
  RAS(0 < newproc && newproc <= ref_mpi_n(ref_mpi),
      "newproc negative or larger then nproc");
  ref_malloc_init(vtx, ref_mpi_n(ref_mpi) + 1, PARM_INT, 0);
  ref_malloc_init(send_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(recv_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  /* use ref_part machinery to set the first global index on the subset parts */
  for (proc = 0; proc < newproc; proc++) {
    vtx[proc] = ref_part_first(ntotal, newproc, proc);
  }
  vtx[newproc] = ntotal;
  /* fill vtx with no global on remaining (unused) parts */
  for (proc = newproc + 1; proc <= ref_mpi_n(ref_mpi); proc++) {
    vtx[proc] = vtx[newproc];
  }
  /* new and old vertex count for my rank */
  nold = (REF_INT)(vtxdist[1 + ref_mpi_rank(ref_mpi)] -
                   vtxdist[ref_mpi_rank(ref_mpi)]);
  nnew = (REF_INT)(vtx[1 + ref_mpi_rank(ref_mpi)] - vtx[ref_mpi_rank(ref_mpi)]);
  ref_malloc_init(part, nnew, PARM_INT, REF_EMPTY);
  ref_malloc_init(xadj, nnew + 1, PARM_INT, 0);
  ref_malloc_init(deg, nold, PARM_INT, 0);
  ref_malloc_init(newdeg, nnew, PARM_INT, 0);
  for (i = 0; i < nold; i++) {
    deg[i] = (REF_INT)(xadjdist[i + 1] - xadjdist[i]);
  }
  for (proc = 0; proc < ref_mpi_n(ref_mpi); proc++) {
    n0 = MAX(vtx[proc], vtxdist[ref_mpi_rank(ref_mpi)]);
    n1 = MIN(vtx[proc + 1], vtxdist[ref_mpi_rank(ref_mpi) + 1]);
    send_size[proc] = (REF_INT)MAX(0, n1 - n0);
  }
  for (proc = 0; proc < ref_mpi_n(ref_mpi); proc++) {
    n0 = MAX(vtx[ref_mpi_rank(ref_mpi)], vtxdist[proc]);
    n1 = MIN(vtx[ref_mpi_rank(ref_mpi) + 1], vtxdist[proc + 1]);
    recv_size[proc] = (REF_INT)MAX(0, n1 - n0);
  }
  RSS(ref_mpi_alltoallv(ref_mpi, deg, send_size, newdeg, recv_size, 1,
                        REF_INT_TYPE),
      "alltoallv degree");
  xadj[0] = 0;
  for (i = 0; i < nnew; i++) {
    xadj[i + 1] = xadj[i] + newdeg[i];
  }
  ref_free(newdeg);
  ref_free(deg);
  for (proc = 0; proc < ref_mpi_n(ref_mpi); proc++) {
    n0 = MAX(vtx[proc], vtxdist[ref_mpi_rank(ref_mpi)]);
    n1 = MIN(vtx[proc + 1], vtxdist[ref_mpi_rank(ref_mpi) + 1]);
    nsend = (REF_INT)MAX(0, n1 - n0);
    send_size[proc] = 0;
    if (0 < nsend) {
      first = (REF_INT)(n0 - vtxdist[ref_mpi_rank(ref_mpi)]);
      send_size[proc] = (REF_INT)(xadjdist[first + nsend] - xadjdist[first]);
    }
  }
  RSS(ref_mpi_alltoall(ref_mpi, send_size, recv_size, REF_INT_TYPE),
      "alltoall sizes");
  nrecv = 0;
  for (proc = 0; proc < ref_mpi_n(ref_mpi); proc++) {
    nrecv += recv_size[proc];
  }
  REIS(nrecv, xadj[nnew], "verify total recv sizes");
  ref_malloc_init(adjncy, nrecv, PARM_INT, 0);
  ref_malloc_init(adjwgt, nrecv, PARM_INT, 0);
  RSS(ref_mpi_alltoallv(ref_mpi, adjncydist, send_size, adjncy, recv_size, 1,
                        REF_INT_TYPE),
      "alltoallv adjncy");
  RSS(ref_mpi_alltoallv(ref_mpi, adjwgtdist, send_size, adjwgt, recv_size, 1,
                        REF_INT_TYPE),
      "alltoallv adjwgt");
  ref_mpi_stopwatch_stop(ref_mpi, "parmetis subset");

  /* split comm and call parmetis */
  RSS(ref_mpi_front_comm(ref_mpi, &split_mpi, newproc), "split comm");
  if (ref_mpi_rank(ref_mpi) < newproc) {
    RSS(ref_migrate_parmetis_wrapper(split_mpi, vtx, xadj, adjncy, adjwgt,
                                     part),
        "parmetis wrapper");
  }
  RSS(ref_mpi_join_comm(split_mpi), "join comm");
  RSS(ref_mpi_free(split_mpi), "free split comm");
  ref_free(adjwgt);
  ref_free(adjncy);
  ref_free(xadj);
  ref_mpi_stopwatch_stop(ref_mpi, "parmetis part");

  /* return part to partdist */
  for (proc = 0; proc < ref_mpi_n(ref_mpi); proc++) {
    n0 = MAX(vtx[ref_mpi_rank(ref_mpi)], vtxdist[proc]);
    n1 = MIN(vtx[ref_mpi_rank(ref_mpi) + 1], vtxdist[proc + 1]);
    send_size[proc] = (REF_INT)MAX(0, n1 - n0);
  }
  for (proc = 0; proc < ref_mpi_n(ref_mpi); proc++) {
    n0 = MAX(vtx[proc], vtxdist[ref_mpi_rank(ref_mpi)]);
    n1 = MIN(vtx[proc + 1], vtxdist[ref_mpi_rank(ref_mpi) + 1]);
    recv_size[proc] = (REF_INT)MAX(0, n1 - n0);
  }

  RSS(ref_mpi_alltoallv(ref_mpi, part, send_size, partdist, recv_size, 1,
                        REF_INT_TYPE),
      "alltoallv adjwgt");
  ref_mpi_stopwatch_stop(ref_mpi, "subset part");

  ref_free(part);
  ref_free(recv_size);
  ref_free(send_size);
  ref_free(vtx);
  return REF_SUCCESS;
}
REF_STATUS ref_migrate_parmetis_part(REF_GRID ref_grid, REF_INT *node_part) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_MIGRATE ref_migrate;
  PARM_INT *vtxdist;
  PARM_INT *xadj, *adjncy, *adjwgt;
  PARM_INT *part;

  REF_GLOB *implied, shift;
  REF_INT node, n, proc, *partition_size, degree;
  REF_INT item, ref;
  REF_INT newpart;

  RSS(ref_node_synchronize_globals(ref_node), "sync global nodes");
  RSS(ref_node_collect_ghost_age(ref_node), "collect ghost age");

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  RSS(ref_migrate_create(&ref_migrate, ref_grid), "create migrate");

  /* skip agglomeration stuff */

  n = 0;
  each_ref_migrate_node(ref_migrate, node) { n++; }

  ref_malloc(partition_size, ref_mpi_n(ref_mpi), REF_INT);
  RSS(ref_mpi_allgather(ref_mpi, &n, partition_size, REF_INT_TYPE),
      "gather size of each part");

  ref_malloc(vtxdist, ref_mpi_n(ref_mpi) + 1, PARM_INT);
  ref_malloc_init(implied, ref_migrate_max(ref_migrate), REF_GLOB, REF_EMPTY);
  ref_malloc(xadj, n + 1, PARM_INT);
  ref_malloc_init(part, n, PARM_INT, ref_mpi_rank(ref_mpi));

  vtxdist[0] = 0;
  each_ref_mpi_part(ref_mpi, proc) {
    vtxdist[proc + 1] = vtxdist[proc] + partition_size[proc];
  }

  shift = vtxdist[ref_mpi_rank(ref_mpi)];
  n = 0;
  xadj[0] = 0;
  each_ref_migrate_node(ref_migrate, node) {
    implied[node] = shift + (REF_GLOB)n;
    RSS(ref_adj_degree(ref_migrate_conn(ref_migrate), node, &degree), "deg");
    RAS(0 < degree, "hanging node island, zero degree");
    xadj[n + 1] = xadj[n] + degree;
    n++;
  }
  RSS(ref_node_ghost_glob(ref_node, implied, 1), "implied ghosts");

  ref_malloc(adjncy, xadj[n], PARM_INT);
  ref_malloc(adjwgt, xadj[n], PARM_INT);

  n = 0;
  each_ref_migrate_node(ref_migrate, node) {
    degree = 0;
    each_ref_adj_node_item_with_ref(ref_migrate_conn(ref_migrate), node, item,
                                    ref) {
      adjncy[xadj[n] + degree] = (PARM_INT)implied[ref];
      adjwgt[xadj[n] + degree] = ref_migrate_age(ref_migrate, node) +
                                 ref_migrate_age(ref_migrate, ref) + 1;
      degree++;
    }
    n++;
  }

  ref_mpi_stopwatch_stop(ref_mpi, "parmetis graph");

  newpart = ref_mpi_n(ref_mpi);
  if (ref_node_n_global(ref_node) < 100000) newpart = 1;

  if (1 == newpart) {
    RSS(ref_migrate_metis_wrapper(ref_mpi, vtxdist, xadj, adjncy, adjwgt, part),
        "metis wrapper");
    ref_mpi_stopwatch_stop(ref_mpi, "metis part");
  } else {
    RSS(ref_migrate_parmetis_subset(ref_mpi, newpart, vtxdist, xadj, adjncy,
                                    adjwgt, part),
        "subset");
  }

  n = 0;
  each_ref_migrate_node(ref_migrate, node) {
    node_part[node] = (REF_INT)part[n];
    n++;
  }

  /* skip agglomeration stuff */

  ref_free(adjwgt);
  ref_free(adjncy);
  ref_free(part);
  ref_free(xadj);
  ref_free(implied);
  ref_free(vtxdist);
  ref_free(partition_size);

  RSS(ref_migrate_report_load_balance(ref_grid, node_part), "report bal");

  RSS(ref_migrate_free(ref_migrate), "free migrate");

  return REF_SUCCESS;
}
#endif

static REF_STATUS ref_migrate_new_part(REF_GRID ref_grid, REF_INT *new_part) {
  if (!ref_mpi_para(ref_grid_mpi(ref_grid))) {
    RSS(ref_migrate_single_part(ref_grid, new_part), "single by nproc");
    return REF_SUCCESS;
  }

  switch (ref_grid_partitioner(ref_grid)) {
    case REF_MIGRATE_SINGLE:
      RSS(ref_migrate_single_part(ref_grid, new_part), "single by method");
      break;
    case REF_MIGRATE_NATIVE_RCB:
      RSS(ref_migrate_native_rcb_part(ref_grid, new_part), "single by method");
      break;
    case REF_MIGRATE_ZOLTAN_GRAPH:
    case REF_MIGRATE_ZOLTAN_RCB:
#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
      RSS(ref_migrate_zoltan_part(ref_grid, new_part), "zoltan part");
      break;
#endif
    case REF_MIGRATE_PARMETIS:
#if defined(HAVE_PARMETIS) && defined(HAVE_MPI)
      RSS(ref_migrate_parmetis_part(ref_grid, new_part), "parmetis part");
      break;
#endif
    case REF_MIGRATE_RECOMMENDED:
#if defined(HAVE_PARMETIS) && defined(HAVE_MPI)
      RSS(ref_migrate_parmetis_part(ref_grid, new_part), "parmetis part");
      break;
#endif
#if !defined(HAVE_PARMETIS) && defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
      RSS(ref_migrate_zoltan_part(ref_grid, new_part), "zoltan part");
      break;
#endif
#if !defined(HAVE_PARMETIS) && !defined(HAVE_ZOLTAN)
      RSS(ref_migrate_native_rcb_part(ref_grid, new_part), "single by method");
      break;
#endif
    default:
      if (ref_grid_once(ref_grid))
        printf(
            "requested partioner method %d"
            " is not recognized or configured\n",
            (int)ref_grid_partitioner(ref_grid));
      RSS(REF_IMPLEMENT, "ref_migrate_method");
      break;
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_migrate_shufflin_node(REF_NODE ref_node) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_GLOB *a_global, *b_global;
  REF_INT part, node;
  REF_INT *a_next;
  REF_DBL *a_real, *b_real;
  REF_DBL *a_aux, *b_aux;
  REF_INT local;
  REF_INT i;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_node_valid_node(ref_node, node) {
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) {
      if (ref_node_part(ref_node, node) < 0 ||
          ref_node_part(ref_node, node) >= ref_mpi_n(ref_mpi)) {
        printf("id %d node %d global " REF_GLOB_FMT " part %d",
               ref_mpi_rank(ref_mpi), node, ref_node_global(ref_node, node),
               ref_node_part(ref_node, node));
        THROW("part out of range");
      }
      a_size[ref_node_part(ref_node, node)]++;
    }
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_total += a_size[part];
  ref_malloc(a_global, a_total, REF_GLOB);
  ref_malloc(a_real, REF_NODE_REAL_PER * a_total, REF_DBL);
  a_aux = NULL;
  if (ref_node_naux(ref_node) > 0)
    ref_malloc(a_aux, ref_node_naux(ref_node) * a_total, REF_DBL);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_total += b_size[part];
  ref_malloc(b_global, b_total, REF_GLOB);
  ref_malloc(b_real, REF_NODE_REAL_PER * b_total, REF_DBL);
  b_aux = NULL;
  if (ref_node_naux(ref_node) > 0)
    ref_malloc(b_aux, ref_node_naux(ref_node) * b_total, REF_DBL);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) a_next[part] =
      a_next[part - 1] + a_size[part - 1];

  each_ref_node_valid_node(ref_node, node) {
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) {
      part = ref_node_part(ref_node, node);
      a_global[a_next[part]] = ref_node_global(ref_node, node);
      for (i = 0; i < REF_NODE_REAL_PER; i++)
        a_real[i + REF_NODE_REAL_PER * a_next[part]] =
            ref_node_real(ref_node, i, node);
      for (i = 0; i < ref_node_naux(ref_node); i++)
        a_aux[i + ref_node_naux(ref_node) * a_next[part]] =
            ref_node_aux(ref_node, i, node);
      a_next[part]++;
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_global, a_size, b_global, b_size, 1,
                        REF_GLOB_TYPE),
      "alltoallv global");

  RSS(ref_mpi_alltoallv(ref_mpi, a_real, a_size, b_real, b_size,
                        REF_NODE_REAL_PER, REF_DBL_TYPE),
      "alltoallv real");

  if (ref_node_naux(ref_node) > 0)
    RSS(ref_mpi_alltoallv(ref_mpi, a_aux, a_size, b_aux, b_size,
                          ref_node_naux(ref_node), REF_DBL_TYPE),
        "alltoallv aux");

  RSS(ref_node_add_many(ref_node, b_total, b_global), "add many");

  for (node = 0; node < b_total; node++) {
    RSS(ref_node_local(ref_node, b_global[node], &local), "local");
    for (i = 0; i < REF_NODE_REAL_PER; i++)
      ref_node_real(ref_node, i, local) = b_real[i + REF_NODE_REAL_PER * node];
    for (i = 0; i < ref_node_naux(ref_node); i++)
      ref_node_aux(ref_node, i, local) =
          b_aux[i + ref_node_naux(ref_node) * node];
    ref_node_part(ref_node, local) = ref_mpi_rank(ref_mpi);
  }

  ref_free(a_next);
  ref_free(b_aux);
  ref_free(b_real);
  ref_free(b_global);
  ref_free(a_aux);
  ref_free(a_real);
  ref_free(a_global);
  ref_free(b_size);
  ref_free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_shufflin_cell(REF_NODE ref_node, REF_CELL ref_cell) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT all_parts[REF_CELL_MAX_SIZE_PER];
  REF_INT nunique;
  REF_INT unique_parts[REF_CELL_MAX_SIZE_PER];
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT part, node, cell, i;
  REF_INT *a_next;
  REF_GLOB *a_c2n, *b_c2n;
  REF_INT *a_parts, *b_parts;
  REF_BOOL need_to_keep;

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      all_parts[node] = ref_node_part(ref_node, nodes[node]);
    }
    RSS(ref_sort_unique_int(ref_cell_node_per(ref_cell), all_parts, &nunique,
                            unique_parts),
        "unique");
    for (node = 0; node < nunique; node++) {
      part = unique_parts[node];
      if (ref_mpi_rank(ref_mpi) != part) a_size[part]++;
    }
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_total += a_size[part];
  ref_malloc(a_c2n, ref_cell_size_per(ref_cell) * a_total, REF_GLOB);
  ref_malloc(a_parts, ref_cell_size_per(ref_cell) * a_total, REF_INT);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_total += b_size[part];
  ref_malloc(b_c2n, ref_cell_size_per(ref_cell) * b_total, REF_GLOB);
  ref_malloc(b_parts, ref_cell_size_per(ref_cell) * b_total, REF_INT);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) a_next[part] =
      a_next[part - 1] + a_size[part - 1];

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      all_parts[node] = ref_node_part(ref_node, nodes[node]);
    }
    RSS(ref_sort_unique_int(ref_cell_node_per(ref_cell), all_parts, &nunique,
                            unique_parts),
        "unique");
    for (node = 0; node < nunique; node++) {
      part = unique_parts[node];
      if (ref_mpi_rank(ref_mpi) != part) {
        for (i = 0; i < ref_cell_node_per(ref_cell); i++) {
          a_c2n[i + ref_cell_size_per(ref_cell) * a_next[part]] =
              ref_node_global(ref_node, nodes[i]);
          a_parts[i + ref_cell_size_per(ref_cell) * a_next[part]] =
              ref_node_part(ref_node, nodes[i]);
        }
        if (ref_cell_last_node_is_an_id(ref_cell)) {
          a_c2n[ref_cell_node_per(ref_cell) +
                ref_cell_size_per(ref_cell) * a_next[part]] =
              (REF_GLOB)nodes[ref_cell_node_per(ref_cell)];
          a_parts[ref_cell_node_per(ref_cell) +
                  ref_cell_size_per(ref_cell) * a_next[part]] = REF_EMPTY;
        }
        a_next[part]++;
      }
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_c2n, a_size, b_c2n, b_size,
                        ref_cell_size_per(ref_cell), REF_GLOB_TYPE),
      "alltoallv c2n");
  RSS(ref_mpi_alltoallv(ref_mpi, a_parts, a_size, b_parts, b_size,
                        ref_cell_size_per(ref_cell), REF_INT_TYPE),
      "alltoallv parts");

  RSS(ref_cell_add_many_global(ref_cell, ref_node, b_total, b_c2n, b_parts,
                               ref_mpi_rank(ref_mpi)),
      "g");

  free(a_next);
  free(b_parts);
  free(b_c2n);
  free(a_parts);
  free(a_c2n);
  free(b_size);
  free(a_size);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    need_to_keep = REF_FALSE;
    for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
      need_to_keep = (need_to_keep || (ref_mpi_rank(ref_mpi) ==
                                       ref_node_part(ref_node, nodes[node])));
    }
    if (!need_to_keep) RSS(ref_cell_remove(ref_cell, cell), "remove");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_migrate_shufflin_geom(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_ADJ ref_adj = ref_geom_adj(ref_geom);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT part, node;
  REF_INT *a_next;
  REF_GLOB *a_int, *b_int;
  REF_DBL *a_real, *b_real;
  REF_INT i, degree, item, geom;
  REF_GLOB global;
  REF_INT descr[REF_GEOM_DESCR_SIZE];

  if (!ref_mpi_para(ref_mpi)) return REF_SUCCESS;

  ref_malloc_init(a_size, ref_mpi_n(ref_mpi), REF_INT, 0);
  ref_malloc_init(b_size, ref_mpi_n(ref_mpi), REF_INT, 0);

  each_ref_node_valid_node(ref_node, node) {
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) {
      RSS(ref_adj_degree(ref_adj, node, &degree), "adj deg");
      a_size[ref_node_part(ref_node, node)] += degree;
    }
  }

  RSS(ref_mpi_alltoall(ref_mpi, a_size, b_size, REF_INT_TYPE),
      "alltoall sizes");

  a_total = 0;
  each_ref_mpi_part(ref_mpi, part) a_total += a_size[part];
  ref_malloc(a_int, REF_GEOM_DESCR_SIZE * a_total, REF_GLOB);
  ref_malloc(a_real, 2 * a_total, REF_DBL);

  b_total = 0;
  each_ref_mpi_part(ref_mpi, part) b_total += b_size[part];
  ref_malloc(b_int, REF_GEOM_DESCR_SIZE * b_total, REF_GLOB);
  ref_malloc(b_real, 2 * b_total, REF_DBL);

  ref_malloc(a_next, ref_mpi_n(ref_mpi), REF_INT);
  a_next[0] = 0;
  each_ref_mpi_worker(ref_mpi, part) {
    a_next[part] = a_next[part - 1] + a_size[part - 1];
  }

  each_ref_node_valid_node(ref_node, node) {
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) {
      each_ref_adj_node_item_with_ref(ref_adj, node, item, geom) {
        part = ref_node_part(ref_node, node);
        each_ref_descr(ref_geom, i) {
          a_int[i + REF_GEOM_DESCR_SIZE * a_next[part]] =
              (REF_GLOB)ref_geom_descr(ref_geom, i, geom);
        }
        a_int[REF_GEOM_DESCR_NODE + REF_GEOM_DESCR_SIZE * a_next[part]] =
            ref_node_global(ref_node, ref_geom_node(ref_geom, geom));
        a_real[0 + 2 * a_next[part]] = ref_geom_param(ref_geom, 0, geom);
        a_real[1 + 2 * a_next[part]] = ref_geom_param(ref_geom, 1, geom);
        a_next[part]++;
      }
    }
  }

  RSS(ref_mpi_alltoallv(ref_mpi, a_int, a_size, b_int, b_size,
                        REF_GEOM_DESCR_SIZE, REF_GLOB_TYPE),
      "alltoallv geom int");
  RSS(ref_mpi_alltoallv(ref_mpi, a_real, a_size, b_real, b_size, 2,
                        REF_DBL_TYPE),
      "alltoallv geom real");

  for (geom = 0; geom < b_total; geom++) {
    each_ref_descr(ref_geom, i) {
      descr[i] = (REF_INT)b_int[i + REF_GEOM_DESCR_SIZE * geom];
    }
    global = b_int[REF_GEOM_DESCR_NODE + REF_GEOM_DESCR_SIZE * geom];
    RSS(ref_node_local(ref_node, global, &node), "g2l");
    descr[REF_GEOM_DESCR_NODE] = node;
    RSS(ref_geom_add_with_descr(ref_geom, descr, &(b_real[2 * geom])),
        "geom add");
  }

  free(a_next);
  free(b_real);
  free(b_int);
  free(a_real);
  free(a_int);
  free(b_size);
  free(a_size);

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_shufflin(REF_GRID ref_grid) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT group, node;
  REF_BOOL need_to_keep;

  if (!ref_mpi_para(ref_grid_mpi(ref_grid))) return REF_SUCCESS;

  RSS(ref_node_synchronize_globals(ref_node), "sync global nodes");

  RSS(ref_migrate_shufflin_node(ref_node), "send out nodes");
  RSS(ref_migrate_shufflin_geom(ref_grid), "geom");

  each_ref_grid_all_ref_cell(ref_grid, group, ref_cell) {
    RSS(ref_migrate_shufflin_cell(ref_node, ref_cell), "cell");
  }

  each_ref_node_valid_node(ref_node, node) {
    if (ref_mpi_rank(ref_mpi) != ref_node_part(ref_node, node)) {
      need_to_keep = REF_FALSE;
      each_ref_grid_2d_3d_ref_cell(ref_grid, group, ref_cell) {
        need_to_keep =
            (need_to_keep || !ref_adj_empty(ref_cell_adj(ref_cell), node));
      }
      if (!need_to_keep) {
        RSS(ref_node_remove_without_global_invalidates_sorted(ref_node, node),
            "remove");
        RSS(ref_geom_remove_all(ref_grid_geom(ref_grid), node), "rm geom");
      }
    }
  }
  RSS(ref_node_rebuild_sorted_global(ref_node), "rebuild");

  RSS(ref_node_ghost_real(ref_node), "ghost real");
  RSS(ref_geom_ghost(ref_grid_geom(ref_grid), ref_node), "ghost geom");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "shuffle");

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_to_balance(REF_GRID ref_grid) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_INT *node_part;

  RSS(ref_node_synchronize_globals(ref_node), "sync global nodes");
  RSS(ref_node_collect_ghost_age(ref_node), "collect ghost age");

  ref_malloc_init(node_part, ref_node_max(ref_node), REF_INT, REF_EMPTY);

  RSS(ref_migrate_new_part(ref_grid, node_part), "new part");

  RSS(ref_node_ghost_int(ref_node, node_part, 1), "ghost part");

  if (NULL != ref_grid_interp(ref_grid)) {
    RSS(ref_interp_from_part(ref_grid_interp(ref_grid), node_part),
        "from part");
  } else {
    for (node = 0; node < ref_node_max(ref_node); node++)
      ref_node_part(ref_node, node) = node_part[node];

    RSS(ref_migrate_shufflin(ref_grid), "shufflin");
  }
  ref_free(node_part);

  return REF_SUCCESS;
}

static REF_ULONG ref_migrate_split_morton(REF_ULONG a) {
  REF_ULONG x = a & 0x1fffff; /* we only look at the first 21 bits */
  x = (x | x << 32) & 0x1f00000000ffff;
  /* shift left 32 bits, OR with self, and
     00011111000000000000000000000000000000001111111111111111 */
  x = (x | x << 16) & 0x1f0000ff0000ff;
  /* shift left 32 bits, OR with self, and
     00011111000000000000000011111111000000000000000011111111 */
  x = (x | x << 8) & 0x100f00f00f00f00f;
  /* shift left 32 bits, OR with self, and
     0001000000001111000000001111000000001111000000001111000000000000 */
  x = (x | x << 4) & 0x10c30c30c30c30c3;
  /* shift left 32 bits, OR with self, and
     0001000011000011000011000011000011000011000011000011000100000000 */
  x = (x | x << 2) & 0x1249249249249249;
  return x;
}
REF_ULONG ref_migrate_morton_id(REF_UINT x, REF_UINT y, REF_UINT z) {
  REF_ULONG answer = 0;
  answer |= ref_migrate_split_morton(x) | ref_migrate_split_morton(y) << 1 |
            ref_migrate_split_morton(z) << 2;
  return answer;
}

REF_STATUS ref_migrate_split_dir(REF_MPI ref_mpi, REF_INT n, REF_DBL *xyz,
                                 REF_INT *dir) {
  REF_DBL mins[3], maxes[3], temp;
  REF_INT i, j;
  *dir = 0;
  if (n == 0) {
    return REF_SUCCESS;
  }
  for (j = 0; j < 3; j++) {
    mins[j] = xyz[j];
    maxes[j] = xyz[j];
  }
  for (i = 1; i < n; i++) {
    for (j = 0; j < 3; j++) {
      mins[j] = MIN(mins[j], xyz[j + 3 * i]);
      maxes[j] = MAX(maxes[j], xyz[j + 3 * i]);
    }
  }
  for (j = 0; j < 3; j++) {
    temp = mins[j];
    RSS(ref_mpi_min(ref_mpi, &temp, &(mins[j]), REF_DBL_TYPE), "min");
    RSS(ref_mpi_bcast(ref_mpi, &(mins[j]), 1, REF_DBL_TYPE), "bcast");
    temp = maxes[j];
    RSS(ref_mpi_max(ref_mpi, &temp, &(maxes[j]), REF_DBL_TYPE), "max");
    RSS(ref_mpi_bcast(ref_mpi, &(maxes[j]), 1, REF_DBL_TYPE), "bcast");
  }
  if ((maxes[1] - mins[1]) >= (maxes[0] - mins[0]) &&
      (maxes[1] - mins[1]) >= (maxes[2] - mins[2]))
    *dir = 1;
  if ((maxes[2] - mins[2]) >= (maxes[0] - mins[0]) &&
      (maxes[2] - mins[2]) >= (maxes[1] - mins[1]))
    *dir = 2;

  return REF_SUCCESS;
}

REF_STATUS ref_migrate_split_ratio(REF_INT number_of_partitions,
                                   REF_DBL *ratio) {
  REF_INT half = number_of_partitions / 2;
  if (ref_math_divisible((REF_DBL)half, (REF_DBL)number_of_partitions)) {
    *ratio = (REF_DBL)half / (REF_DBL)number_of_partitions;
  } else {
    *ratio = 0;
    return REF_DIV_ZERO;
  }
  return REF_SUCCESS;
}
