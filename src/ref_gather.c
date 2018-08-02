
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
#include <string.h>

#include "ref_export.h"
#include "ref_gather.h"

#include "ref_endian.h"
#include "ref_malloc.h"
#include "ref_mpi.h"
#include "ref_sort.h"

REF_STATUS ref_gather_create(REF_GATHER *ref_gather_ptr) {
  REF_GATHER ref_gather;

  ref_malloc(*ref_gather_ptr, 1, REF_GATHER_STRUCT);

  ref_gather = *ref_gather_ptr;

  ref_gather->recording = REF_FALSE;
  ref_gather->file = (FILE *)NULL;
  ref_gather->time = 0.0;

  return REF_SUCCESS;
}

REF_STATUS ref_gather_free(REF_GATHER ref_gather) {
  if (NULL != (void *)(ref_gather->file)) fclose(ref_gather->file);
  ref_free(ref_gather);

  return REF_SUCCESS;
}

REF_STATUS ref_gather_tec_movie_record_button(REF_GATHER ref_gather,
                                              REF_BOOL on_or_off) {
  ref_gather->recording = on_or_off;
  return REF_SUCCESS;
}

static REF_STATUS ref_gather_ncell_below_quality(REF_NODE ref_node,
                                                 REF_CELL ref_cell,
                                                 REF_DBL min_quality,
                                                 REF_INT *ncell) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ncell_local;
  REF_DBL quality;

  ncell_local = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    if (ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, nodes[0]) &&
        quality < min_quality)
      ncell_local++;
  }

  RSS(ref_mpi_sum(ref_mpi, &ncell_local, ncell, 1, REF_INT_TYPE), "sum");
  RSS(ref_mpi_bcast(ref_mpi, ncell, 1, REF_INT_TYPE), "bcast");

  return REF_SUCCESS;
}

static REF_STATUS ref_gather_node_tec_part(REF_NODE ref_node, REF_INT nnode,
                                           REF_INT *l2c, REF_DBL *scalar,
                                           FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT chunk;
  REF_DBL *local_xyzm, *xyzm;
  REF_INT nnode_written, first, n, i;
  REF_INT global, local;
  REF_STATUS status;
  REF_INT dim = 6;
  REF_INT *sorted_local, *sorted_cellnode, *pack, total_cellnode, position;

  total_cellnode = 0;
  for (i = 0; i < ref_node_max(ref_node); i++) {
    if (REF_EMPTY != l2c[i] && ref_node_owned(ref_node, i)) {
      total_cellnode++;
    }
  }

  ref_malloc(sorted_local, total_cellnode, REF_INT);
  ref_malloc(sorted_cellnode, total_cellnode, REF_INT);
  ref_malloc(pack, total_cellnode, REF_INT);

  total_cellnode = 0;
  for (i = 0; i < ref_node_max(ref_node); i++) {
    if (REF_EMPTY != l2c[i] && ref_node_owned(ref_node, i)) {
      sorted_cellnode[total_cellnode] = l2c[i];
      pack[total_cellnode] = i;
      total_cellnode++;
    }
  }
  RSS(ref_sort_heap_int(total_cellnode, sorted_cellnode, sorted_local), "sort");
  for (i = 0; i < total_cellnode; i++) {
    sorted_local[i] = pack[sorted_local[i]];
    sorted_cellnode[i] = l2c[sorted_local[i]];
  }
  ref_free(pack);

  chunk = nnode / ref_mpi_n(ref_mpi) + 1;
  chunk = MAX(chunk, 100000);

  ref_malloc(local_xyzm, dim * chunk, REF_DBL);
  ref_malloc(xyzm, dim * chunk, REF_DBL);

  nnode_written = 0;
  while (nnode_written < nnode) {
    first = nnode_written;
    n = MIN(chunk, nnode - nnode_written);

    nnode_written += n;

    for (i = 0; i < dim * chunk; i++) local_xyzm[i] = 0.0;

    for (i = 0; i < n; i++) {
      global = first + i;
      status =
          ref_sort_search(total_cellnode, sorted_cellnode, global, &position);
      RXS(status, REF_NOT_FOUND, "node local failed");
      if (REF_SUCCESS == status) {
        local = sorted_local[position];
        local_xyzm[0 + dim * i] = ref_node_xyz(ref_node, 0, local);
        local_xyzm[1 + dim * i] = ref_node_xyz(ref_node, 1, local);
        local_xyzm[2 + dim * i] = ref_node_xyz(ref_node, 2, local);
        local_xyzm[3 + dim * i] = (REF_DBL)ref_node_part(ref_node, local);
        if (NULL == scalar) {
          local_xyzm[4 + dim * i] = (REF_DBL)ref_node_age(ref_node, local);
        } else {
          local_xyzm[4 + dim * i] = scalar[local];
        }
        local_xyzm[5 + dim * i] = 1.0;
      } else {
        local_xyzm[0 + dim * i] = 0.0;
        local_xyzm[1 + dim * i] = 0.0;
        local_xyzm[2 + dim * i] = 0.0;
        local_xyzm[3 + dim * i] = 0.0;
        local_xyzm[4 + dim * i] = 0.0;
        local_xyzm[5 + dim * i] = 0.0;
      }
    }

    for (i = 0; i < n; i++)
      if ((ABS(local_xyzm[5 + dim * i] - 1.0) > 0.1) &&
          (ABS(local_xyzm[5 + dim * i] - 0.0) > 0.1)) {
        printf("%s: %d: %s: before sum %d %f\n", __FILE__, __LINE__, __func__,
               first + i, local_xyzm[5 + dim * i]);
      }

    RSS(ref_mpi_sum(ref_mpi, local_xyzm, xyzm, dim * n, REF_DBL_TYPE), "sum");

    if (ref_mpi_once(ref_mpi))
      for (i = 0; i < n; i++) {
        if (ABS(xyzm[5 + dim * i] - 1.0) > 0.1) {
          printf("%s: %d: %s: after sum %d %f\n", __FILE__, __LINE__, __func__,
                 first + i, xyzm[5 + dim * i]);
        }
        if (NULL == scalar) {
          fprintf(file, "%.15e %.15e %.15e %.0f %.0f\n", xyzm[0 + dim * i],
                  xyzm[1 + dim * i], xyzm[2 + dim * i], xyzm[3 + dim * i],
                  xyzm[4 + dim * i]);
        } else {
          fprintf(file, "%.15e %.15e %.15e %.0f %.8f\n", xyzm[0 + dim * i],
                  xyzm[1 + dim * i], xyzm[2 + dim * i], xyzm[3 + dim * i],
                  xyzm[4 + dim * i]);
        }
      }
  }

  ref_free(xyzm);
  ref_free(local_xyzm);
  ref_free(sorted_cellnode);
  ref_free(sorted_local);

  return REF_SUCCESS;
}

static REF_STATUS ref_gather_cell_tec(REF_NODE ref_node, REF_CELL ref_cell,
                                      REF_INT ncell_expected, REF_INT *l2c,
                                      FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT cell, node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per = ref_cell_node_per(ref_cell);
  REF_INT *c2n, ncell;
  REF_INT proc, part;
  REF_INT ncell_actual;

  ncell_actual = 0;

  if (ref_mpi_once(ref_mpi)) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "part");
      if (ref_mpi_rank(ref_mpi) == part) {
        for (node = 0; node < node_per; node++) {
          nodes[node] = l2c[nodes[node]];
          nodes[node]++;
          fprintf(file, " %d", nodes[node]);
        }
        ncell_actual++;
        fprintf(file, "\n");
      }
    }
  }

  if (ref_mpi_once(ref_mpi)) {
    each_ref_mpi_worker(ref_mpi, proc) {
      RSS(ref_mpi_recv(ref_mpi, &ncell, 1, REF_INT_TYPE, proc), "recv ncell");
      ref_malloc(c2n, ncell * node_per, REF_INT);
      RSS(ref_mpi_recv(ref_mpi, c2n, ncell * node_per, REF_INT_TYPE, proc),
          "recv c2n");
      for (cell = 0; cell < ncell; cell++) {
        for (node = 0; node < node_per; node++) {
          c2n[node + node_per * cell]++;
          fprintf(file, " %d", c2n[node + node_per * cell]);
        }
        ncell_actual++;
        fprintf(file, "\n");
      }
      ref_free(c2n);
    }
  } else {
    ncell = 0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "part");
      if (ref_mpi_rank(ref_mpi) == part) ncell++;
    }
    RSS(ref_mpi_send(ref_mpi, &ncell, 1, REF_INT_TYPE, 0), "send ncell");
    ref_malloc(c2n, ncell * node_per, REF_INT);
    ncell = 0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "part");
      if (ref_mpi_rank(ref_mpi) == part) {
        for (node = 0; node < node_per; node++)
          c2n[node + node_per * ncell] = l2c[nodes[node]];
        ncell++;
      }
    }
    RSS(ref_mpi_send(ref_mpi, c2n, ncell * node_per, REF_INT_TYPE, 0),
        "send c2n");

    ref_free(c2n);
  }

  if (ref_mpi_once(ref_mpi)) {
    REIS(ncell_expected, ncell_actual, "cell count mismatch");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_gather_cell_quality_tec(REF_NODE ref_node,
                                              REF_CELL ref_cell,
                                              REF_DBL min_quality, FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT cell, node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per = ref_cell_node_per(ref_cell);
  REF_INT size_per = ref_cell_size_per(ref_cell);
  REF_INT ncell;
  REF_INT *c2n;
  REF_INT proc;
  REF_DBL quality;

  if (ref_mpi_once(ref_mpi)) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
      if (ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, nodes[0]) &&
          quality < min_quality) {
        for (node = 0; node < node_per; node++) {
          nodes[node] = ref_node_global(ref_node, nodes[node]);
          nodes[node]++;
          fprintf(file, " %d", nodes[node]);
        }
        fprintf(file, "\n");
      }
    }
  }

  if (ref_mpi_once(ref_mpi)) {
    each_ref_mpi_worker(ref_mpi, proc) {
      RSS(ref_mpi_recv(ref_mpi, &ncell, 1, REF_INT_TYPE, proc), "recv ncell");
      if (ncell > 0) {
        ref_malloc(c2n, ncell * size_per, REF_INT);
        RSS(ref_mpi_recv(ref_mpi, c2n, ncell * size_per, REF_INT_TYPE, proc),
            "recv c2n");
        for (cell = 0; cell < ncell; cell++) {
          for (node = 0; node < node_per; node++) {
            c2n[node + size_per * cell]++;
            fprintf(file, " %d", c2n[node + size_per * cell]);
          }
          fprintf(file, "\n");
        }
        ref_free(c2n);
      }
    }
  } else {
    ncell = 0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
      if (ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, nodes[0]) &&
          quality < min_quality)
        ncell++;
    }
    RSS(ref_mpi_send(ref_mpi, &ncell, 1, REF_INT_TYPE, 0), "send ncell");
    if (ncell > 0) {
      ref_malloc(c2n, ncell * size_per, REF_INT);
      ncell = 0;
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
        if (ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, nodes[0]) &&
            quality < min_quality) {
          for (node = 0; node < node_per; node++)
            c2n[node + size_per * ncell] =
                ref_node_global(ref_node, nodes[node]);
          for (node = node_per; node < size_per; node++)
            c2n[node + size_per * ncell] = nodes[node];
          ncell++;
        }
      }
      RSS(ref_mpi_send(ref_mpi, c2n, ncell * size_per, REF_INT_TYPE, 0),
          "send c2n");
      ref_free(c2n);
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_gather_plt_movie_frame(REF_GRID ref_grid,
                                      const char *zone_title) {
  REF_GATHER ref_gather = ref_grid_gather(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_INT nnode, ncell, *l2c;

  int one = 1;
  int filetype = 0; /* 0 = FULL,1 = GRID,2 = SOLUTION */
  int ascii[1024];
  int numvar = 5;

  float zonemarker = 299.0;
  int parentzone = -1;
  int strandid = -1;
  double solutiontime = (double)(ref_gather->time);
  int notused = -1;
  int zonetype = 2;    /* 0=ORDERED,1=FELINESEG,2=FETRIANGLE, 3=FEQUADRILATERAL,
                        * 4=FETETRAHEDRON, 5=FEBRICK, 6=FEPOLYGON,7=FEPOLYHEDRON */
  int datapacking = 0; /*0=Block, point does not work.*/
  int varloc = 0;      /*0 = Don't specify, all data is located at nodes*/
  int faceneighbors = 0;
  int numpts;
  int numelements;
  int celldim = 0;
  int aux = 0;
  float eohmarker = 357.0;
  int dataformat = 1; /* 1=Float,2=Double,3=LongInt,4=ShortInt,5=Byte,6=Bit */
  int passive = 0;
  int varsharing = 0;
  int connsharing = -1;

  REF_DBL mindata, maxdata, reduction;
  double var_limit;

  REF_INT i, length, ixyz, node;
  REF_DBL *nodal;
  float nodal_float;
  int c2n;
  REF_INT cell, part, incoming, *conn;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  if (!(ref_gather->recording)) return REF_SUCCESS;

  RSS(ref_node_synchronize_globals(ref_node), "sync");

  RSS(ref_grid_cell_nodes(ref_grid, ref_cell, &nnode, &ncell, &l2c), "compact");
  if (ref_grid_once(ref_grid)) {
    if (NULL == (void *)(ref_gather->file)) {
      ref_gather->file = fopen("ref_gather_movie.plt", "w");
      if (NULL == (void *)(ref_gather->file))
        printf("unable to open ref_gather_movie.plt\n");
      RNS(ref_gather->file, "unable to open file");

      REIS(8, fwrite(&"#!TDV112", sizeof(char), 8, ref_gather->file), "header");
      REIS(1, fwrite(&one, sizeof(int), 1, ref_gather->file), "magic");
      REIS(1, fwrite(&filetype, sizeof(int), 1, ref_gather->file), "filetype");

      ascii[0] = (int)'f';
      ascii[1] = (int)'t';
      ascii[2] = 0;
      REIS(3, fwrite(&ascii, sizeof(int), 3, ref_gather->file), "title");

      REIS(1, fwrite(&numvar, sizeof(int), 1, ref_gather->file), "numvar");
      ascii[0] = (int)'x';
      ascii[1] = 0;
      REIS(2, fwrite(&ascii, sizeof(int), 2, ref_gather->file), "var");
      ascii[0] = (int)'y';
      ascii[1] = 0;
      REIS(2, fwrite(&ascii, sizeof(int), 2, ref_gather->file), "var");
      ascii[0] = (int)'z';
      ascii[1] = 0;
      REIS(2, fwrite(&ascii, sizeof(int), 2, ref_gather->file), "var");
      ascii[0] = (int)'p';
      ascii[1] = 0;
      REIS(2, fwrite(&ascii, sizeof(int), 2, ref_gather->file), "var");
      ascii[0] = (int)'a';
      ascii[1] = 0;
      REIS(2, fwrite(&ascii, sizeof(int), 2, ref_gather->file), "var");
    }

    REIS(1, fwrite(&zonemarker, sizeof(float), 1, ref_gather->file),
         "zonemarker");

    if (NULL == zone_title) {
      ascii[0] = (int)'p';
      ascii[1] = (int)'a';
      ascii[2] = (int)'r';
      ascii[3] = (int)'t';
      ascii[4] = 0;
      length = 5;
    } else {
      length = strlen(zone_title);
      RAS(length < 1023, "title too long");
      for (i = 0; i < length; i++) {
        ascii[i] = (int)zone_title[i];
      }
      ascii[length] = 0;
      length++;
    }
    REIS(length, fwrite(&ascii, sizeof(int), length, ref_gather->file),
         "title");

    REIS(1, fwrite(&parentzone, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&strandid, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&solutiontime, sizeof(double), 1, ref_gather->file),
         "double");
    REIS(1, fwrite(&notused, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&zonetype, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&datapacking, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&varloc, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&faceneighbors, sizeof(int), 1, ref_gather->file), "int");
    numpts = nnode;
    REIS(1, fwrite(&numpts, sizeof(int), 1, ref_gather->file), "int");
    numelements = ncell;
    REIS(1, fwrite(&numelements, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&celldim, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&celldim, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&celldim, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&aux, sizeof(int), 1, ref_gather->file), "int");

    REIS(1, fwrite(&eohmarker, sizeof(float), 1, ref_gather->file),
         "eohmarker");
    REIS(1, fwrite(&zonemarker, sizeof(float), 1, ref_gather->file),
         "zonemarker");

    REIS(1, fwrite(&dataformat, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&dataformat, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&dataformat, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&dataformat, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&dataformat, sizeof(int), 1, ref_gather->file), "int");

    REIS(1, fwrite(&passive, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&varsharing, sizeof(int), 1, ref_gather->file), "int");
    REIS(1, fwrite(&connsharing, sizeof(int), 1, ref_gather->file), "int");
  }

  for (ixyz = 0; ixyz < 3; ixyz++) {
    mindata = 1.0e100;
    maxdata = -1.0e100;
    each_ref_node_valid_node(ref_node, node) {
      if (REF_EMPTY != l2c[node]) {
        mindata = MIN(mindata, ref_node_xyz(ref_node, ixyz, node));
        maxdata = MAX(maxdata, ref_node_xyz(ref_node, ixyz, node));
      }
    }
    RSS(ref_mpi_min(ref_mpi, &mindata, &reduction, REF_DBL_TYPE), "min");
    RSS(ref_mpi_bcast(ref_mpi, &reduction, 1, REF_DBL_TYPE), "bcast");
    var_limit = reduction;
    if (ref_grid_once(ref_grid))
      REIS(1, fwrite(&var_limit, sizeof(double), 1, ref_gather->file), "write");
    RSS(ref_mpi_max(ref_mpi, &maxdata, &reduction, REF_DBL_TYPE), "max");
    RSS(ref_mpi_bcast(ref_mpi, &reduction, 1, REF_DBL_TYPE), "bcast");
    var_limit = reduction;
    if (ref_grid_once(ref_grid))
      REIS(1, fwrite(&var_limit, sizeof(double), 1, ref_gather->file), "write");
  }

  mindata = 1.0e100;
  maxdata = -1.0e100;
  each_ref_node_valid_node(ref_node, node) {
    if (REF_EMPTY != l2c[node]) {
      mindata = MIN(mindata, (REF_DBL)ref_node_part(ref_node, node));
      maxdata = MAX(maxdata, (REF_DBL)ref_node_part(ref_node, node));
    }
  }
  RSS(ref_mpi_min(ref_mpi, &mindata, &reduction, REF_DBL_TYPE), "min");
  RSS(ref_mpi_bcast(ref_mpi, &reduction, 1, REF_DBL_TYPE), "bcast");
  var_limit = reduction;
  if (ref_grid_once(ref_grid))
    REIS(1, fwrite(&var_limit, sizeof(double), 1, ref_gather->file), "write");
  RSS(ref_mpi_max(ref_mpi, &maxdata, &reduction, REF_DBL_TYPE), "max");
  RSS(ref_mpi_bcast(ref_mpi, &reduction, 1, REF_DBL_TYPE), "bcast");
  var_limit = reduction;
  if (ref_grid_once(ref_grid))
    REIS(1, fwrite(&var_limit, sizeof(double), 1, ref_gather->file), "write");

  mindata = 1.0e100;
  maxdata = -1.0e100;
  each_ref_node_valid_node(ref_node, node) {
    if (REF_EMPTY != l2c[node]) {
      mindata = MIN(mindata, (REF_DBL)ref_node_age(ref_node, node));
      maxdata = MAX(maxdata, (REF_DBL)ref_node_age(ref_node, node));
    }
  }
  RSS(ref_mpi_min(ref_mpi, &mindata, &reduction, REF_DBL_TYPE), "min");
  RSS(ref_mpi_bcast(ref_mpi, &reduction, 1, REF_DBL_TYPE), "bcast");
  var_limit = reduction;
  if (ref_grid_once(ref_grid))
    REIS(1, fwrite(&var_limit, sizeof(double), 1, ref_gather->file), "write");
  RSS(ref_mpi_max(ref_mpi, &maxdata, &reduction, REF_DBL_TYPE), "max");
  RSS(ref_mpi_bcast(ref_mpi, &reduction, 1, REF_DBL_TYPE), "bcast");
  var_limit = reduction;
  if (ref_grid_once(ref_grid))
    REIS(1, fwrite(&var_limit, sizeof(double), 1, ref_gather->file), "write");

  ref_malloc_init(nodal, nnode * 6, REF_DBL, 0.0);
  each_ref_node_valid_node(ref_node, node) {
    if (REF_EMPTY != l2c[node] && ref_node_owned(ref_node, node)) {
      nodal[0 + 6 * l2c[node]] = ref_node_xyz(ref_node, 0, node);
      nodal[1 + 6 * l2c[node]] = ref_node_xyz(ref_node, 1, node);
      nodal[2 + 6 * l2c[node]] = ref_node_xyz(ref_node, 2, node);
      nodal[3 + 6 * l2c[node]] = (REF_DBL)ref_node_part(ref_node, node);
      nodal[4 + 6 * l2c[node]] = (REF_DBL)ref_node_age(ref_node, node);
      nodal[5 + 6 * l2c[node]] = 1.0;
    }
  }
  RSS(ref_mpi_allsum(ref_mpi, nodal, nnode * 6, REF_DBL_TYPE), "allsum");
  for (node = 0; node < nnode; node++) {
    RWDS(1.0, nodal[5 + 6 * node], -1.0, "node contribution");
  }
  if (ref_grid_once(ref_grid)) {
    for (i = 0; i < 5; i++) {
      for (node = 0; node < nnode; node++) {
        nodal_float = (float)nodal[i + 6 * node];
        REIS(1, fwrite(&nodal_float, sizeof(float), 1, ref_gather->file),
             "nodal");
      }
    }
  }
  ref_free(nodal);
  if (ref_grid_once(ref_grid)) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "part");
      if (ref_mpi_rank(ref_mpi) == part) {
        for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
          c2n = l2c[nodes[node]];
          REIS(1, fwrite(&c2n, sizeof(int), 1, ref_gather->file), "nodal");
        }
      }
    }
    each_ref_mpi_worker(ref_mpi, part) {
      RSS(ref_mpi_recv(ref_mpi, &incoming, 1, REF_INT_TYPE, part), "in");
      if (incoming > 0) {
        ref_malloc(conn, ref_cell_node_per(ref_cell) * incoming, REF_INT);
        RSS(ref_mpi_recv(ref_mpi, &conn, ref_cell_node_per(ref_cell) * incoming,
                         REF_INT_TYPE, part),
            "conn");
        for (node = 0; node < ref_cell_node_per(ref_cell) * incoming; node++) {
          c2n = conn[node];
          REIS(1, fwrite(&c2n, sizeof(int), 1, ref_gather->file), "nodal");
        }
        ref_free(conn);
      }
    }
  } else {
    incoming = 0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "part");
      if (ref_mpi_rank(ref_mpi) == part) {
        incoming++;
      }
    }
    RSS(ref_mpi_send(ref_mpi, &incoming, 1, REF_INT_TYPE, 0), "in");
    if (incoming > 0) {
      ref_malloc(conn, ref_cell_node_per(ref_cell) * incoming, REF_INT);
      each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
        RSS(ref_cell_part(ref_cell, ref_node, cell, &part), "part");
        if (ref_mpi_rank(ref_mpi) == part) {
          for (node = 0; node < ref_cell_node_per(ref_cell); node++) {
            conn[node + ref_cell_node_per(ref_cell) * incoming] =
                l2c[nodes[node]];
          }
        }
        RSS(ref_mpi_send(ref_mpi, &conn, ref_cell_node_per(ref_cell) * incoming,
                         REF_INT_TYPE, 0),
            "conn");
        ref_free(conn);
      }
    }
  }
  ref_free(l2c);

  return REF_SUCCESS;
}

REF_STATUS ref_gather_tec_movie_frame(REF_GRID ref_grid,
                                      const char *zone_title) {
  REF_GATHER ref_gather = ref_grid_gather(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT nnode, ncell, *l2c;
  REF_DBL *norm_dev, dot;

  if (!(ref_gather->recording)) return REF_SUCCESS;

  norm_dev = NULL;
  if (ref_geom_model_loaded(ref_geom)) {
    ref_malloc_init(norm_dev, ref_node_max(ref_node), REF_DBL, 2.0);
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      RSS(ref_geom_tri_norm_deviation(ref_grid, nodes, &dot), "norm dev");
      each_ref_cell_cell_node(ref_cell, cell_node) {
        norm_dev[nodes[cell_node]] = MIN(norm_dev[nodes[cell_node]], dot);
      }
    }
  }

  RSS(ref_node_synchronize_globals(ref_node), "sync");

  RSS(ref_grid_cell_nodes(ref_grid, ref_cell, &nnode, &ncell, &l2c), "l2c");

  if (ref_grid_once(ref_grid)) {
    if (NULL == (void *)(ref_gather->file)) {
      ref_gather->file = fopen("ref_gather_movie.tec", "w");
      if (NULL == (void *)(ref_gather->file))
        printf("unable to open ref_gather_movie.tec\n");
      RNS(ref_gather->file, "unable to open file");

      fprintf(ref_gather->file, "title=\"tecplot refine partion file\"\n");
      fprintf(ref_gather->file, "variables = \"x\" \"y\" \"z\" \"p\" \"a\"\n");
    }
    if (NULL == zone_title) {
      fprintf(ref_gather->file,
              "zone t=\"part\", nodes=%d, elements=%d, datapacking=%s, "
              "zonetype=%s, solutiontime=%f\n",
              nnode, ncell, "point", "fetriangle", ref_gather->time);
    } else {
      fprintf(ref_gather->file,
              "zone t=\"%s\", nodes=%d, elements=%d, datapacking=%s, "
              "zonetype=%s, solutiontime=%f\n",
              zone_title, nnode, ncell, "point", "fetriangle",
              ref_gather->time);
    }
  }

  RSS(ref_gather_node_tec_part(ref_node, nnode, l2c, norm_dev,
                               ref_gather->file),
      "nodes");
  RSS(ref_gather_cell_tec(ref_node, ref_cell, ncell, l2c, ref_gather->file),
      "t");

  ref_free(norm_dev);
  norm_dev = NULL;

  if (REF_FALSE) {
    REF_INT ntet;
    REF_DBL min_quality = 0.10;
    RSS(ref_gather_ncell_below_quality(ref_node, ref_grid_tet(ref_grid),
                                       min_quality, &ntet),
        "ntri");

    if (ref_grid_once(ref_grid)) {
      if (NULL == zone_title) {
        fprintf(ref_gather->file,
                "zone t=\"qpart\", nodes=%d, elements=%d, datapacking=%s, "
                "zonetype=%s, solutiontime=%f\n",
                nnode, MAX(1, ntet), "point", "fetetrahedron",
                ref_gather->time);
      } else {
        fprintf(ref_gather->file,
                "zone t=\"q%s\", nodes=%d, elements=%d, datapacking=%s, "
                "zonetype=%s, solutiontime=%f\n",
                zone_title, nnode, MAX(1, ntet), "point", "fetetrahedron",
                ref_gather->time);
      }
    }
    RSS(ref_gather_node_tec_part(ref_node, nnode, l2c, NULL, ref_gather->file),
        "nodes");
    if (0 == ntet) {
      if (ref_grid_once(ref_grid)) fprintf(ref_gather->file, " 1 1 1 1\n");
    } else {
      RSS(ref_gather_cell_quality_tec(ref_node, ref_grid_tet(ref_grid),
                                      min_quality, ref_gather->file),
          "qtet");
    }
  }

  if (ref_grid_once(ref_grid)) {
    REIS(0, fflush(ref_gather->file), "gather movie fflush");
    (ref_gather->time) += 1.0;
  }

  ref_free(l2c);

  return REF_SUCCESS;
}

REF_STATUS ref_gather_tec_part(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT nnode, ncell, *l2c;

  RSS(ref_node_synchronize_globals(ref_node), "sync");

  RSS(ref_grid_cell_nodes(ref_grid, ref_cell, &nnode, &ncell, &l2c), "l2c");

  file = NULL;
  if (ref_grid_once(ref_grid)) {
    file = fopen(filename, "w");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    fprintf(file, "title=\"tecplot refine partion file\"\n");
    fprintf(file, "variables = \"x\" \"y\" \"z\" \"p\" \"a\"\n");
    fprintf(
        file,
        "zone t=\"part\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
        nnode, ncell, "point", "fetriangle");
  }

  RSS(ref_gather_node_tec_part(ref_node, nnode, l2c, NULL, file), "nodes");
  RSS(ref_gather_cell_tec(ref_node, ref_cell, ncell, l2c, file), "nodes");

  if (ref_grid_once(ref_grid)) fclose(file);

  ref_free(l2c);

  return REF_SUCCESS;
}

static REF_STATUS ref_gather_node(REF_NODE ref_node, REF_BOOL swap_endian,
                                  REF_BOOL has_id, FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT chunk;
  REF_DBL *local_xyzm, *xyzm;
  REF_DBL swapped_dbl;
  REF_INT nnode_written, first, n, i;
  REF_INT global, local;
  REF_INT id = REF_EXPORT_MESHB_VERTEX_ID;
  REF_STATUS status;

  chunk = ref_node_n_global(ref_node) / ref_mpi_n(ref_mpi) + 1;

  ref_malloc(local_xyzm, 4 * chunk, REF_DBL);
  ref_malloc(xyzm, 4 * chunk, REF_DBL);

  nnode_written = 0;
  while (nnode_written < ref_node_n_global(ref_node)) {
    first = nnode_written;
    n = MIN(chunk, ref_node_n_global(ref_node) - nnode_written);

    nnode_written += n;

    for (i = 0; i < 4 * chunk; i++) local_xyzm[i] = 0.0;

    for (i = 0; i < n; i++) {
      global = first + i;
      status = ref_node_local(ref_node, global, &local);
      RXS(status, REF_NOT_FOUND, "node local failed");
      if (REF_SUCCESS == status &&
          ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, local)) {
        local_xyzm[0 + 4 * i] = ref_node_xyz(ref_node, 0, local);
        local_xyzm[1 + 4 * i] = ref_node_xyz(ref_node, 1, local);
        local_xyzm[2 + 4 * i] = ref_node_xyz(ref_node, 2, local);
        local_xyzm[3 + 4 * i] = 1.0;
      } else {
        local_xyzm[0 + 4 * i] = 0.0;
        local_xyzm[1 + 4 * i] = 0.0;
        local_xyzm[2 + 4 * i] = 0.0;
        local_xyzm[3 + 4 * i] = 0.0;
      }
    }

    RSS(ref_mpi_sum(ref_mpi, local_xyzm, xyzm, 4 * n, REF_DBL_TYPE), "sum");

    if (ref_mpi_once(ref_mpi))
      for (i = 0; i < n; i++) {
        if (ABS(xyzm[3 + 4 * i] - 1.0) > 0.1) {
          printf("error gather node %d %f\n", first + i, xyzm[3 + 4 * i]);
        }
        swapped_dbl = xyzm[0 + 4 * i];
        if (swap_endian) SWAP_DBL(swapped_dbl);
        REIS(1, fwrite(&swapped_dbl, sizeof(REF_DBL), 1, file), "x");
        swapped_dbl = xyzm[1 + 4 * i];
        if (swap_endian) SWAP_DBL(swapped_dbl);
        REIS(1, fwrite(&swapped_dbl, sizeof(REF_DBL), 1, file), "y");
        swapped_dbl = xyzm[2 + 4 * i];
        if (swap_endian) SWAP_DBL(swapped_dbl);
        REIS(1, fwrite(&swapped_dbl, sizeof(REF_DBL), 1, file), "z");
        if (has_id) REIS(1, fwrite(&id, sizeof(REF_INT), 1, file), "id");
      }
  }

  ref_free(xyzm);
  ref_free(local_xyzm);

  return REF_SUCCESS;
}

static REF_STATUS ref_gather_node_metric(REF_NODE ref_node, FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT chunk;
  REF_DBL *local_xyzm, *xyzm;
  REF_INT nnode_written, first, n, i, im;
  REF_INT global, local;
  REF_STATUS status;

  chunk = ref_node_n_global(ref_node) / ref_mpi_n(ref_mpi) + 1;

  ref_malloc(local_xyzm, 7 * chunk, REF_DBL);
  ref_malloc(xyzm, 7 * chunk, REF_DBL);

  nnode_written = 0;
  while (nnode_written < ref_node_n_global(ref_node)) {
    first = nnode_written;
    n = MIN(chunk, ref_node_n_global(ref_node) - nnode_written);

    nnode_written += n;

    for (i = 0; i < 7 * chunk; i++) local_xyzm[i] = 0.0;

    for (i = 0; i < n; i++) {
      global = first + i;
      status = ref_node_local(ref_node, global, &local);
      RXS(status, REF_NOT_FOUND, "node local failed");
      if (REF_SUCCESS == status &&
          ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, local)) {
        for (im = 0; im < 6; im++)
          local_xyzm[im + 7 * i] = ref_node_metric(ref_node, im, local);
        local_xyzm[6 + 7 * i] = 1.0;
      } else {
        for (im = 0; im < 7; im++) local_xyzm[im + 7 * i] = 0.0;
      }
    }

    RSS(ref_mpi_sum(ref_mpi, local_xyzm, xyzm, 7 * n, REF_DBL_TYPE), "sum");

    if (ref_mpi_once(ref_mpi))
      for (i = 0; i < n; i++) {
        if (ABS(xyzm[6 + 7 * i] - 1.0) > 0.1) {
          printf("error gather node %d %f\n", first + i, xyzm[6 + 7 * i]);
        }
        fprintf(file, "%.15e %.15e %.15e %.15e %.15e %.15e \n", xyzm[0 + 7 * i],
                xyzm[1 + 7 * i], xyzm[2 + 7 * i], xyzm[3 + 7 * i],
                xyzm[4 + 7 * i], xyzm[5 + 7 * i]);
      }
  }

  ref_free(xyzm);
  ref_free(local_xyzm);

  return REF_SUCCESS;
}

static REF_STATUS ref_gather_node_metric_solb(REF_NODE ref_node, FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT chunk;
  REF_DBL *local_xyzm, *xyzm;
  REF_INT nnode_written, first, n, i, im;
  REF_INT global, local;
  REF_STATUS status;
  int next_position, keyword_code;
  int code, version, dim;

  if (ref_mpi_once(ref_mpi)) {
    code = 1;
    REIS(1, fwrite(&code, sizeof(int), 1, file), "code");
    version = 2;
    REIS(1, fwrite(&version, sizeof(int), 1, file), "version");
    next_position = 4 + 4 + 4 + ftell(file);
    keyword_code = 3;
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "dim code");
    REIS(1, fwrite(&next_position, sizeof(next_position), 1, file), "next pos");
    dim = 3;
    REIS(1, fwrite(&dim, sizeof(int), 1, file), "dim");
  }

  if (ref_mpi_once(ref_mpi)) {
    next_position =
        4 + 4 + 4 + 4 + 4 + ref_node_n_global(ref_node) * (6 * 8) + ftell(file);
    keyword_code = 62;
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "vertex version code");
    REIS(1, fwrite(&next_position, sizeof(next_position), 1, file), "next pos");
    REIS(1, fwrite(&(ref_node_n_global(ref_node)), sizeof(int), 1, file),
         "nnode");
    keyword_code = 1; /* one solution at node */
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "n solutions");
    keyword_code = 3; /* solution type 3, metric */
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "metric solution");
  }

  chunk = ref_node_n_global(ref_node) / ref_mpi_n(ref_mpi) + 1;

  ref_malloc(local_xyzm, 7 * chunk, REF_DBL);
  ref_malloc(xyzm, 7 * chunk, REF_DBL);

  nnode_written = 0;
  while (nnode_written < ref_node_n_global(ref_node)) {
    first = nnode_written;
    n = MIN(chunk, ref_node_n_global(ref_node) - nnode_written);

    nnode_written += n;

    for (i = 0; i < 7 * chunk; i++) local_xyzm[i] = 0.0;

    for (i = 0; i < n; i++) {
      global = first + i;
      status = ref_node_local(ref_node, global, &local);
      RXS(status, REF_NOT_FOUND, "node local failed");
      if (REF_SUCCESS == status &&
          ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, local)) {
        for (im = 0; im < 6; im++)
          local_xyzm[im + 7 * i] = ref_node_metric(ref_node, im, local);
        local_xyzm[6 + 7 * i] = 1.0;
      } else {
        for (im = 0; im < 7; im++) local_xyzm[im + 7 * i] = 0.0;
      }
    }

    RSS(ref_mpi_sum(ref_mpi, local_xyzm, xyzm, 7 * n, REF_DBL_TYPE), "sum");

    if (ref_mpi_once(ref_mpi))
      for (i = 0; i < n; i++) {
        if (ABS(xyzm[6 + 7 * i] - 1.0) > 0.1) {
          printf("error gather node %d %f\n", first + i, xyzm[6 + 7 * i]);
        }
        REIS(1, fwrite(&(xyzm[0 + 7 * i]), sizeof(REF_DBL), 1, file), "m11");
        REIS(1, fwrite(&(xyzm[1 + 7 * i]), sizeof(REF_DBL), 1, file), "m12");
        /* transposed 3,2 */
        REIS(1, fwrite(&(xyzm[3 + 7 * i]), sizeof(REF_DBL), 1, file), "m22");
        REIS(1, fwrite(&(xyzm[2 + 7 * i]), sizeof(REF_DBL), 1, file), "m13");
        REIS(1, fwrite(&(xyzm[4 + 7 * i]), sizeof(REF_DBL), 1, file), "m23");
        REIS(1, fwrite(&(xyzm[5 + 7 * i]), sizeof(REF_DBL), 1, file), "m33");
      }
  }

  ref_free(xyzm);
  ref_free(local_xyzm);

  if (ref_mpi_once(ref_mpi))
    REIS(next_position, ftell(file), "solb metric record len inconsistent");

  if (ref_mpi_once(ref_mpi)) { /* End */
    keyword_code = 54;
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "end kw");
    next_position = 0;
    REIS(1, fwrite(&next_position, sizeof(next_position), 1, file), "next pos");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_gather_node_scalar_solb(REF_NODE ref_node, REF_INT ldim,
                                              REF_DBL *scalar, FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT chunk;
  REF_DBL *local_xyzm, *xyzm;
  REF_INT nnode_written, first, n, i, im;
  REF_INT global, local;
  REF_STATUS status;
  int next_position, keyword_code;
  int code, version, dim;

  if (ref_mpi_once(ref_mpi)) {
    code = 1;
    REIS(1, fwrite(&code, sizeof(int), 1, file), "code");
    version = 2;
    REIS(1, fwrite(&version, sizeof(int), 1, file), "version");
    next_position = 4 + 4 + 4 + ftell(file);
    keyword_code = 3;
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "dim code");
    REIS(1, fwrite(&next_position, sizeof(next_position), 1, file), "next pos");
    dim = 3;
    REIS(1, fwrite(&dim, sizeof(int), 1, file), "dim");
  }

  if (ref_mpi_once(ref_mpi)) {
    next_position = 4 + 4 + 4 + 4 + 4 +
                    ref_node_n_global(ref_node) * (ldim * 8) + ftell(file);
    keyword_code = 62;
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "vertex version code");
    REIS(1, fwrite(&next_position, sizeof(next_position), 1, file), "next pos");
    REIS(1, fwrite(&(ref_node_n_global(ref_node)), sizeof(int), 1, file),
         "nnode");
    keyword_code = ldim; /* one solution at node */
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "n solutions");
    keyword_code = 1; /* solution type 1, scalar */
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "metric solution");
  }

  chunk = ref_node_n_global(ref_node) / ref_mpi_n(ref_mpi) + 1;

  ref_malloc(local_xyzm, (ldim + 1) * chunk, REF_DBL);
  ref_malloc(xyzm, (ldim + 1) * chunk, REF_DBL);

  nnode_written = 0;
  while (nnode_written < ref_node_n_global(ref_node)) {
    first = nnode_written;
    n = MIN(chunk, ref_node_n_global(ref_node) - nnode_written);

    nnode_written += n;

    for (i = 0; i < (ldim + 1) * chunk; i++) local_xyzm[i] = 0.0;

    for (i = 0; i < n; i++) {
      global = first + i;
      status = ref_node_local(ref_node, global, &local);
      RXS(status, REF_NOT_FOUND, "node local failed");
      if (REF_SUCCESS == status &&
          ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, local)) {
        for (im = 0; im < ldim; im++)
          local_xyzm[im + (ldim + 1) * i] = scalar[im + ldim * local];
        local_xyzm[ldim + (ldim + 1) * i] = 1.0;
      } else {
        for (im = 0; im < (ldim + 1); im++)
          local_xyzm[im + (ldim + 1) * i] = 0.0;
      }
    }

    RSS(ref_mpi_sum(ref_mpi, local_xyzm, xyzm, (ldim + 1) * n, REF_DBL_TYPE),
        "sum");

    if (ref_mpi_once(ref_mpi))
      for (i = 0; i < n; i++) {
        if (ABS(xyzm[ldim + (ldim + 1) * i] - 1.0) > 0.1) {
          printf("error gather node %d %f\n", first + i,
                 xyzm[ldim + (ldim + 1) * i]);
        }
        for (im = 0; im < ldim; im++) {
          REIS(1,
               fwrite(&(xyzm[im + (ldim + 1) * i]), sizeof(REF_DBL), 1, file),
               "s");
        }
      }
  }

  ref_free(xyzm);
  ref_free(local_xyzm);

  if (ref_mpi_once(ref_mpi))
    REIS(next_position, ftell(file), "solb metric record len inconsistent");

  if (ref_mpi_once(ref_mpi)) { /* End */
    keyword_code = 54;
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "end kw");
    next_position = 0;
    REIS(1, fwrite(&next_position, sizeof(next_position), 1, file), "next pos");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_gather_cell(REF_NODE ref_node, REF_CELL ref_cell,
                                  REF_BOOL faceid_insted_of_c2n,
                                  REF_BOOL always_id, REF_BOOL swap_endian,
                                  REF_BOOL select_faceid, REF_INT faceid,
                                  FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT cell, node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node_per = ref_cell_node_per(ref_cell);
  REF_INT size_per = ref_cell_size_per(ref_cell);
  REF_INT ncell;
  REF_INT *c2n;
  REF_INT proc;

  if (ref_mpi_once(ref_mpi)) {
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, nodes[0]) &&
          (!select_faceid || nodes[ref_cell_node_per(ref_cell)] == faceid)) {
        if (faceid_insted_of_c2n) {
          node = node_per;
          if (swap_endian) SWAP_INT(nodes[node]);
          REIS(1, fwrite(&(nodes[node]), sizeof(REF_INT), 1, file), "cel node");
        } else {
          for (node = 0; node < node_per; node++) {
            nodes[node] = ref_node_global(ref_node, nodes[node]);
            nodes[node]++;
            if (swap_endian) SWAP_INT(nodes[node]);
            REIS(1, fwrite(&(nodes[node]), sizeof(REF_INT), 1, file),
                 "cel node");
          }
          if (always_id) {
            if (ref_cell_last_node_is_an_id(ref_cell)) {
              node = node_per;
              if (swap_endian) SWAP_INT(nodes[node]);
              REIS(1, fwrite(&(nodes[node]), sizeof(REF_INT), 1, file),
                   "cel node");
            } else {
              node = REF_EXPORT_MESHB_3D_ID;
              if (swap_endian) SWAP_INT(node);
              REIS(1, fwrite(&(node), sizeof(REF_INT), 1, file), "cel node");
            }
          }
        }
      }
    }
  }

  if (ref_mpi_once(ref_mpi)) {
    each_ref_mpi_worker(ref_mpi, proc) {
      RSS(ref_mpi_recv(ref_mpi, &ncell, 1, REF_INT_TYPE, proc), "recv ncell");
      if (ncell > 0) {
        ref_malloc(c2n, ncell * size_per, REF_INT);
        RSS(ref_mpi_recv(ref_mpi, c2n, ncell * size_per, REF_INT_TYPE, proc),
            "recv c2n");
        for (cell = 0; cell < ncell; cell++)
          if (faceid_insted_of_c2n) {
            node = node_per;
            if (swap_endian) SWAP_INT(c2n[node + size_per * cell]);
            REIS(1,
                 fwrite(&(c2n[node + size_per * cell]), sizeof(REF_INT), 1,
                        file),
                 "cell");
          } else {
            for (node = 0; node < node_per; node++) {
              c2n[node + size_per * cell]++;
              if (swap_endian) SWAP_INT(c2n[node + size_per * cell]);
              REIS(1,
                   fwrite(&(c2n[node + size_per * cell]), sizeof(REF_INT), 1,
                          file),
                   "cell");
            }
            if (always_id) {
              if (ref_cell_last_node_is_an_id(ref_cell)) {
                node = node_per;
                if (swap_endian) SWAP_INT(c2n[node + size_per * cell]);
                REIS(1,
                     fwrite(&(c2n[node + size_per * cell]), sizeof(REF_INT), 1,
                            file),
                     "cel node");
              } else {
                node = REF_EXPORT_MESHB_3D_ID;
                if (swap_endian) SWAP_INT(node);
                REIS(1, fwrite(&(node), sizeof(REF_INT), 1, file), "cel node");
              }
            }
          }
        ref_free(c2n);
      }
    }
  } else {
    ncell = 0;
    each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
      if (ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, nodes[0]) &&
          (!select_faceid || nodes[ref_cell_node_per(ref_cell)] == faceid))
        ncell++;
      RSS(ref_mpi_send(ref_mpi, &ncell, 1, REF_INT_TYPE, 0), "send ncell");
      if (ncell > 0) {
        ref_malloc(c2n, ncell * size_per, REF_INT);
        ncell = 0;
        each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
          if (ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, nodes[0]) &&
              (!select_faceid ||
               nodes[ref_cell_node_per(ref_cell)] == faceid)) {
            for (node = 0; node < node_per; node++)
              c2n[node + size_per * ncell] =
                  ref_node_global(ref_node, nodes[node]);
            for (node = node_per; node < size_per; node++)
              c2n[node + size_per * ncell] = nodes[node];
            ncell++;
          }
        }
        RSS(ref_mpi_send(ref_mpi, c2n, ncell * size_per, REF_INT_TYPE, 0),
            "send c2n");
        ref_free(c2n);
      }
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_gather_geom(REF_NODE ref_node, REF_GEOM ref_geom,
                                  REF_INT type, FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT geom, node, id, i;
  REF_INT ngeom;
  REF_INT *node_id;
  REF_DBL *param;
  REF_INT proc;
  double filler = 0.0;

  if (ref_mpi_once(ref_mpi)) {
    each_ref_geom_of(ref_geom, type, geom) {
      if (ref_mpi_rank(ref_mpi) !=
          ref_node_part(ref_node, ref_geom_node(ref_geom, geom)))
        continue;
      node = ref_node_global(ref_node, ref_geom_node(ref_geom, geom)) + 1;
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
  }

  if (ref_mpi_once(ref_mpi)) {
    each_ref_mpi_worker(ref_mpi, proc) {
      RSS(ref_mpi_recv(ref_mpi, &ngeom, 1, REF_INT_TYPE, proc), "recv ngeom");
      if (ngeom > 0) {
        ref_malloc(node_id, 2 * ngeom, REF_INT);
        ref_malloc(param, 2 * ngeom, REF_DBL);
        RSS(ref_mpi_recv(ref_mpi, node_id, 2 * ngeom, REF_INT_TYPE, proc),
            "recv node_id");
        RSS(ref_mpi_recv(ref_mpi, param, 2 * ngeom, REF_DBL_TYPE, proc),
            "recv param");
        for (geom = 0; geom < ngeom; geom++) {
          node = node_id[0 + 2 * geom] + 1;
          id = node_id[1 + 2 * geom];
          REIS(1, fwrite(&(node), sizeof(int), 1, file), "node");
          REIS(1, fwrite(&(id), sizeof(int), 1, file), "id");
          for (i = 0; i < type; i++)
            REIS(1, fwrite(&(param[i + 2 * geom]), sizeof(double), 1, file),
                 "id");
          if (0 < type)
            REIS(1, fwrite(&(filler), sizeof(double), 1, file), "id");
        }
        ref_free(param);
        ref_free(node_id);
      }
    }
  } else {
    ngeom = 0;
    each_ref_geom_of(ref_geom, type, geom) {
      if (ref_mpi_rank(ref_mpi) !=
          ref_node_part(ref_node, ref_geom_node(ref_geom, geom)))
        continue;
      ngeom++;
    }
    RSS(ref_mpi_send(ref_mpi, &ngeom, 1, REF_INT_TYPE, 0), "send ngeom");
    if (ngeom > 0) {
      ref_malloc(node_id, 2 * ngeom, REF_INT);
      ref_malloc_init(param, 2 * ngeom, REF_DBL, 0.0); /* prevent uninit */
      ngeom = 0;
      each_ref_geom_of(ref_geom, type, geom) {
        if (ref_mpi_rank(ref_mpi) !=
            ref_node_part(ref_node, ref_geom_node(ref_geom, geom)))
          continue;
        node_id[0 + 2 * ngeom] =
            ref_node_global(ref_node, ref_geom_node(ref_geom, geom));
        node_id[1 + 2 * ngeom] = ref_geom_id(ref_geom, geom);
        for (i = 0; i < type; i++)
          param[i + 2 * ngeom] = ref_geom_param(ref_geom, i, geom);
        ngeom++;
      }
      RSS(ref_mpi_send(ref_mpi, node_id, 2 * ngeom, REF_INT_TYPE, 0),
          "send node_id");
      RSS(ref_mpi_send(ref_mpi, param, 2 * ngeom, REF_DBL_TYPE, 0),
          "send param");
      ref_free(param);
      ref_free(node_id);
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_gather_meshb(REF_GRID ref_grid, const char *filename) {
  REF_BOOL verbose = REF_FALSE;
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT code, version, dim;
  int next_position, keyword_code;
  REF_INT ncell, node_per;
  REF_INT ngeom, type;
  REF_CELL ref_cell;
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_BOOL faceid_insted_of_c2n = REF_FALSE;
  REF_BOOL always_id = REF_TRUE;
  REF_BOOL swap_endian = REF_FALSE;
  REF_BOOL select_faceid = REF_FALSE;
  REF_INT faceid = REF_EMPTY;

  RAS(!ref_grid_twod(ref_grid), "only 3D");

  RSS(ref_node_synchronize_globals(ref_node), "sync");
  file = NULL;
  if (ref_grid_once(ref_grid)) {
    file = fopen(filename, "w");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    code = 1;
    REIS(1, fwrite(&code, sizeof(int), 1, file), "code");
    version = 2;
    REIS(1, fwrite(&version, sizeof(int), 1, file), "version");
    next_position = 4 + 4 + 4 + ftell(file);
    keyword_code = 3;
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "dim code");
    REIS(1, fwrite(&next_position, sizeof(next_position), 1, file), "next pos");
    dim = 3;
    REIS(1, fwrite(&dim, sizeof(int), 1, file), "dim");
  }

  if (ref_grid_once(ref_grid)) {
    next_position =
        4 + 4 + 4 + ref_node_n_global(ref_node) * (3 * 8 + 4) + ftell(file);
    keyword_code = 4;
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "vertex version code");
    REIS(1, fwrite(&next_position, sizeof(next_position), 1, file), "next pos");
    REIS(1, fwrite(&(ref_node_n_global(ref_node)), sizeof(int), 1, file),
         "nnode");
    if (verbose)
      printf("vertex kw %d next %d n %d\n", keyword_code, next_position,
             ref_node_n_global(ref_node));
  }
  RSS(ref_gather_node(ref_node, swap_endian, always_id, file), "nodes");
  if (ref_grid_once(ref_grid))
    REIS(next_position, ftell(file), "vertex inconsistent");

  ref_cell = ref_grid_edg(ref_grid);
  keyword_code = 5;
  RSS(ref_gather_ncell(ref_node, ref_cell, &ncell), "ntet");
  if (verbose)
    printf("%d: edge ncell %d\n", ref_mpi_rank(ref_grid_mpi(ref_grid)), ncell);
  if (ncell > 0) {
    if (ref_grid_once(ref_grid)) {
      node_per = ref_cell_node_per(ref_cell);
      next_position = 4 + 4 + 4 + ncell * (4 * (node_per + 1)) + ftell(file);
      REIS(1, fwrite(&keyword_code, sizeof(int), 1, file),
           "vertex version code");
      REIS(1, fwrite(&next_position, sizeof(next_position), 1, file),
           "next pos");
      REIS(1, fwrite(&(ncell), sizeof(int), 1, file), "nnode");
      if (verbose)
        printf("elem kw %d next %d n %d\n", keyword_code, next_position, ncell);
    }
    RSS(ref_gather_cell(ref_node, ref_cell, faceid_insted_of_c2n, always_id,
                        swap_endian, select_faceid, faceid, file),
        "nodes");
    if (ref_grid_once(ref_grid))
      REIS(next_position, ftell(file), "cell inconsistent");
  }

  ref_cell = ref_grid_tri(ref_grid);
  keyword_code = 6;
  RSS(ref_gather_ncell(ref_node, ref_cell, &ncell), "ntet");
  if (ncell > 0) {
    if (ref_grid_once(ref_grid)) {
      node_per = ref_cell_node_per(ref_cell);
      next_position = 4 + 4 + 4 + ncell * (4 * (node_per + 1)) + ftell(file);
      REIS(1, fwrite(&keyword_code, sizeof(int), 1, file),
           "vertex version code");
      REIS(1, fwrite(&next_position, sizeof(next_position), 1, file),
           "next pos");
      REIS(1, fwrite(&(ncell), sizeof(int), 1, file), "nnode");
      if (verbose)
        printf("elem kw %d next %d n %d\n", keyword_code, next_position, ncell);
    }
    RSS(ref_gather_cell(ref_node, ref_cell, faceid_insted_of_c2n, always_id,
                        swap_endian, select_faceid, faceid, file),
        "nodes");
    if (ref_grid_once(ref_grid))
      REIS(next_position, ftell(file), "cell inconsistent");
  }

  ref_cell = ref_grid_tet(ref_grid);
  keyword_code = 8;
  RSS(ref_gather_ncell(ref_node, ref_cell, &ncell), "ntet");
  if (ncell > 0) {
    if (ref_grid_once(ref_grid)) {
      node_per = ref_cell_node_per(ref_cell);
      next_position = 4 + 4 + 4 + ncell * (4 * (node_per + 1)) + ftell(file);
      REIS(1, fwrite(&keyword_code, sizeof(int), 1, file),
           "vertex version code");
      REIS(1, fwrite(&next_position, sizeof(next_position), 1, file),
           "next pos");
      REIS(1, fwrite(&(ncell), sizeof(int), 1, file), "nnode");
      if (verbose)
        printf("elem kw %d next %d n %d\n", keyword_code, next_position, ncell);
    }
    RSS(ref_gather_cell(ref_node, ref_cell, faceid_insted_of_c2n, always_id,
                        swap_endian, select_faceid, faceid, file),
        "nodes");
    if (ref_grid_once(ref_grid))
      REIS(next_position, ftell(file), "cell inconsistent");
  }

  each_ref_type(ref_geom, type) {
    keyword_code = 40 + type; /* GmfVerticesOnGeometricVertices */
    RSS(ref_gather_ngeom(ref_node, ref_geom, type, &ngeom), "ngeom");
    if (ngeom > 0) {
      if (ref_grid_once(ref_grid)) {
        node_per = ref_cell_node_per(ref_cell);
        next_position = 4 + 4 + 4 + ngeom * (4 * 2 + 8 * type) +
                        (0 < type ? 8 * ngeom : 0) + ftell(file);
        REIS(1, fwrite(&keyword_code, sizeof(int), 1, file),
             "vertex version code");
        REIS(1, fwrite(&next_position, sizeof(next_position), 1, file),
             "next pos");
        REIS(1, fwrite(&(ngeom), sizeof(int), 1, file), "nnode");
        if (verbose)
          printf("geom type %d kw %d next %d n %d\n", type, keyword_code,
                 next_position, ngeom);
      }
      RSS(ref_gather_geom(ref_node, ref_geom, type, file), "nodes");
      if (ref_grid_once(ref_grid))
        REIS(next_position, ftell(file), "cell inconsistent");
    }
  }

  if (ref_grid_once(ref_grid) && 0 < ref_geom_cad_data_size(ref_geom)) {
    keyword_code = 126; /* GmfByteFlow */
    next_position = 4 + 4 + 4 + ref_geom_cad_data_size(ref_geom) + ftell(file);
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "keyword");
    REIS(1, fwrite(&next_position, sizeof(next_position), 1, file), "next");
    REIS(1, fwrite(&(ref_geom_cad_data_size(ref_geom)), sizeof(int), 1, file),
         "n");
    REIS(ref_geom_cad_data_size(ref_geom),
         fwrite(ref_geom_cad_data(ref_geom), sizeof(REF_BYTE),
                ref_geom_cad_data_size(ref_geom), file),
         "node");
    REIS(next_position, ftell(file), "cad_model inconsistent");
  }

  if (ref_grid_once(ref_grid)) { /* End */
    keyword_code = 54;           /* GmfEnd 101-47 */
    REIS(1, fwrite(&keyword_code, sizeof(int), 1, file), "vertex version code");
    next_position = 0;
    REIS(1, fwrite(&next_position, sizeof(next_position), 1, file), "next pos");
    if (verbose) printf("end kw %d next %d\n", keyword_code, next_position);
    if (verbose) printf("close %s\n", filename);
    fclose(file);
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_gather_bin_ugrid(REF_GRID ref_grid, const char *filename,
                                       REF_BOOL swap_endian) {
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nnode, ntri, nqua, ntet, npyr, npri, nhex;
  REF_CELL ref_cell;
  REF_INT group;
  REF_INT faceid, min_faceid, max_faceid;
  REF_BOOL always_id = REF_FALSE;
  REF_BOOL faceid_insted_of_c2n, select_faceid;

  RSS(ref_node_synchronize_globals(ref_node), "sync");

  nnode = ref_node_n_global(ref_node);

  RSS(ref_gather_ncell(ref_node, ref_grid_tri(ref_grid), &ntri), "ntri");
  RSS(ref_gather_ncell(ref_node, ref_grid_qua(ref_grid), &nqua), "nqua");

  RSS(ref_gather_ncell(ref_node, ref_grid_tet(ref_grid), &ntet), "ntet");
  RSS(ref_gather_ncell(ref_node, ref_grid_pyr(ref_grid), &npyr), "npyr");
  RSS(ref_gather_ncell(ref_node, ref_grid_pri(ref_grid), &npri), "npri");
  RSS(ref_gather_ncell(ref_node, ref_grid_hex(ref_grid), &nhex), "nhex");

  file = NULL;
  if (ref_grid_once(ref_grid)) {
    file = fopen(filename, "w");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    if (swap_endian) {
      SWAP_INT(nnode);
      SWAP_INT(ntri);
      SWAP_INT(nqua);
      SWAP_INT(ntet);
      SWAP_INT(npyr);
      SWAP_INT(npri);
      SWAP_INT(nhex);
    }

    REIS(1, fwrite(&nnode, sizeof(REF_INT), 1, file), "nnode");

    REIS(1, fwrite(&ntri, sizeof(REF_INT), 1, file), "ntri");
    REIS(1, fwrite(&nqua, sizeof(REF_INT), 1, file), "nqua");

    REIS(1, fwrite(&ntet, sizeof(REF_INT), 1, file), "ntet");
    REIS(1, fwrite(&npyr, sizeof(REF_INT), 1, file), "npyr");
    REIS(1, fwrite(&npri, sizeof(REF_INT), 1, file), "npri");
    REIS(1, fwrite(&nhex, sizeof(REF_INT), 1, file), "nhex");
  }

  RSS(ref_gather_node(ref_node, swap_endian, always_id, file), "nodes");

  RSS(ref_geom_faceid_range(ref_grid, &min_faceid, &max_faceid), "range");

  faceid_insted_of_c2n = REF_FALSE;
  select_faceid = REF_TRUE;
  for (faceid = min_faceid; faceid <= max_faceid; faceid++)
    RSS(ref_gather_cell(ref_node, ref_grid_tri(ref_grid), faceid_insted_of_c2n,
                        always_id, swap_endian, select_faceid, faceid, file),
        "tri c2n");
  for (faceid = min_faceid; faceid <= max_faceid; faceid++)
    RSS(ref_gather_cell(ref_node, ref_grid_qua(ref_grid), faceid_insted_of_c2n,
                        always_id, swap_endian, select_faceid, faceid, file),
        "qua c2n");

  faceid_insted_of_c2n = REF_TRUE;
  for (faceid = min_faceid; faceid <= max_faceid; faceid++)
    RSS(ref_gather_cell(ref_node, ref_grid_tri(ref_grid), faceid_insted_of_c2n,
                        always_id, swap_endian, select_faceid, faceid, file),
        "tri faceid");
  for (faceid = min_faceid; faceid <= max_faceid; faceid++)
    RSS(ref_gather_cell(ref_node, ref_grid_qua(ref_grid), faceid_insted_of_c2n,
                        always_id, swap_endian, select_faceid, faceid, file),
        "qua faceid");

  faceid_insted_of_c2n = REF_FALSE;
  select_faceid = REF_FALSE;
  faceid = REF_EMPTY;
  each_ref_grid_ref_cell(ref_grid, group, ref_cell) {
    RSS(ref_gather_cell(ref_node, ref_cell, faceid_insted_of_c2n, always_id,
                        swap_endian, select_faceid, faceid, file),
        "cell c2n");
  }

  if (ref_grid_once(ref_grid)) fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_gather_by_extension(REF_GRID ref_grid, const char *filename) {
  size_t end_of_string;

  end_of_string = strlen(filename);

  if (strcmp(&filename[end_of_string - 10], ".lb8.ugrid") == 0) {
    RSS(ref_gather_bin_ugrid(ref_grid, filename, REF_FALSE),
        "lb8_ugrid failed");
    return REF_SUCCESS;
  }
  if (strcmp(&filename[end_of_string - 9], ".b8.ugrid") == 0) {
    RSS(ref_gather_bin_ugrid(ref_grid, filename, REF_TRUE), "b8_ugrid failed");
    return REF_SUCCESS;
  }
  if (strcmp(&filename[end_of_string - 6], ".meshb") == 0) {
    RSS(ref_gather_meshb(ref_grid, filename), "meshb failed");
    return REF_SUCCESS;
  }
  printf("%s: %d: %s %s\n", __FILE__, __LINE__,
         "input file name extension unknown", filename);
  RSS(REF_FAILURE, "unknown file extension");
  return REF_FAILURE;
}

REF_STATUS ref_gather_metric(REF_GRID ref_grid, const char *filename) {
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  size_t end_of_string;
  REF_BOOL solb_format = REF_FALSE;

  RSS(ref_node_synchronize_globals(ref_node), "sync");

  file = NULL;
  if (ref_grid_once(ref_grid)) {
    file = fopen(filename, "w");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    end_of_string = strlen(filename);
    if (strcmp(&filename[end_of_string - 5], ".solb") == 0)
      solb_format = REF_TRUE;
  }
  RSS(ref_mpi_all_or(ref_grid_mpi(ref_grid), &solb_format), "bcast");

  if (solb_format) {
    RSS(ref_gather_node_metric_solb(ref_node, file), "nodes");
  } else {
    RSS(ref_gather_node_metric(ref_node, file), "nodes");
  }

  if (ref_grid_once(ref_grid)) fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_gather_scalar(REF_GRID ref_grid, REF_INT ldim, REF_DBL *scalar,
                             const char *filename) {
  FILE *file;
  REF_NODE ref_node = ref_grid_node(ref_grid);

  RSS(ref_node_synchronize_globals(ref_node), "sync");

  file = NULL;
  if (ref_grid_once(ref_grid)) {
    file = fopen(filename, "w");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");
  }

  RSS(ref_gather_node_scalar_solb(ref_node, ldim, scalar, file), "nodes");

  if (ref_grid_once(ref_grid)) fclose(file);

  return REF_SUCCESS;
}

REF_STATUS ref_gather_ncell(REF_NODE ref_node, REF_CELL ref_cell,
                            REF_INT *ncell) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ncell_local;

  ncell_local = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, nodes[0])) {
      ncell_local++;
    }
  }

  RSS(ref_mpi_sum(ref_mpi, &ncell_local, ncell, 1, REF_INT_TYPE), "sum");
  RSS(ref_mpi_bcast(ref_mpi, ncell, 1, REF_INT_TYPE), "bcast");

  return REF_SUCCESS;
}

REF_STATUS ref_gather_ngeom(REF_NODE ref_node, REF_GEOM ref_geom, REF_INT type,
                            REF_INT *ngeom) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT geom, node;
  REF_INT ngeom_local;

  ngeom_local = 0;
  each_ref_geom_of(ref_geom, type, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (ref_mpi_rank(ref_mpi) == ref_node_part(ref_node, node)) ngeom_local++;
  }

  RSS(ref_mpi_sum(ref_mpi, &ngeom_local, ngeom, 1, REF_INT_TYPE), "sum");
  RSS(ref_mpi_bcast(ref_mpi, ngeom, 1, REF_INT_TYPE), "bcast");

  return REF_SUCCESS;
}
