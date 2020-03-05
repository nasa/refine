
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

#include "ref_part.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_dict.h"
#include "ref_endian.h"
#include "ref_export.h"
#include "ref_import.h"
#include "ref_malloc.h"
#include "ref_migrate.h"
#include "ref_mpi.h"

static REF_STATUS ref_part_node(FILE *file, REF_BOOL swap_endian,
                                REF_BOOL has_id, REF_BOOL twod,
                                REF_NODE ref_node, REF_LONG nnode) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT node, new_node;
  REF_INT part;
  REF_INT n, id;
  REF_DBL dbl;
  REF_DBL *xyz;

  RSS(ref_node_initialize_n_global(ref_node, nnode), "init nnodesg");

  if (ref_mpi_once(ref_mpi)) {
    part = 0;
    for (node = 0; node < ref_part_first(nnode, ref_mpi_n(ref_mpi), 1);
         node++) {
      RSS(ref_node_add(ref_node, node, &new_node), "new_node");
      ref_node_part(ref_node, new_node) = ref_mpi_rank(ref_mpi);
      RES(1, fread(&dbl, sizeof(REF_DBL), 1, file), "x");
      if (swap_endian) SWAP_DBL(dbl);
      ref_node_xyz(ref_node, 0, new_node) = dbl;
      RES(1, fread(&dbl, sizeof(REF_DBL), 1, file), "y");
      if (swap_endian) SWAP_DBL(dbl);
      ref_node_xyz(ref_node, 1, new_node) = dbl;
      if (!twod) {
        RES(1, fread(&dbl, sizeof(REF_DBL), 1, file), "z");
        if (swap_endian) SWAP_DBL(dbl);
      } else {
        dbl = 0.0;
      }
      ref_node_xyz(ref_node, 2, new_node) = dbl;
      if (has_id) REIS(1, fread(&(id), sizeof(id), 1, file), "id");
    }
    each_ref_mpi_worker(ref_mpi, part) {
      n = (REF_INT)(ref_part_first(nnode, ref_mpi_n(ref_mpi), part + 1) -
                    ref_part_first(nnode, ref_mpi_n(ref_mpi), part));
      RSS(ref_mpi_send(ref_mpi, &n, 1, REF_INT_TYPE, part), "send");
      if (n > 0) {
        ref_malloc(xyz, 3 * n, REF_DBL);
        for (node = 0; node < n; node++) {
          RES(1, fread(&dbl, sizeof(REF_DBL), 1, file), "x");
          if (swap_endian) SWAP_DBL(dbl);
          xyz[0 + 3 * node] = dbl;
          RES(1, fread(&dbl, sizeof(REF_DBL), 1, file), "y");
          if (swap_endian) SWAP_DBL(dbl);
          xyz[1 + 3 * node] = dbl;
          if (!twod) {
            RES(1, fread(&dbl, sizeof(REF_DBL), 1, file), "z");
            if (swap_endian) SWAP_DBL(dbl);
          } else {
            dbl = 0.0;
          }
          xyz[2 + 3 * node] = dbl;
          if (has_id) REIS(1, fread(&(id), sizeof(id), 1, file), "id");
        }
        RSS(ref_mpi_send(ref_mpi, xyz, 3 * n, REF_DBL_TYPE, part), "send");
        free(xyz);
      }
    }
  } else {
    RSS(ref_mpi_recv(ref_mpi, &n, 1, REF_INT_TYPE, 0), "recv");
    if (n > 0) {
      ref_malloc(xyz, 3 * n, REF_DBL);
      RSS(ref_mpi_recv(ref_mpi, xyz, 3 * n, REF_DBL_TYPE, 0), "recv");
      for (node = 0; node < n; node++) {
        RSS(ref_node_add(ref_node,
                         node + ref_part_first(nnode, ref_mpi_n(ref_mpi),
                                               ref_mpi_rank(ref_mpi)),
                         &new_node),
            "new_node");
        ref_node_part(ref_node, new_node) = ref_mpi_rank(ref_mpi);
        ref_node_xyz(ref_node, 0, new_node) = xyz[0 + 3 * node];
        ref_node_xyz(ref_node, 1, new_node) = xyz[1 + 3 * node];
        ref_node_xyz(ref_node, 2, new_node) = xyz[2 + 3 * node];
      }
      free(xyz);
    }
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_part_meshb_geom(REF_GEOM ref_geom, REF_INT ngeom,
                                      REF_INT type, REF_NODE ref_node,
                                      REF_INT nnode, FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT end_of_message = REF_EMPTY;
  REF_INT chunk;
  REF_INT *sent_node;
  REF_INT *sent_id;
  REF_DBL *sent_param;
  REF_INT *read_node;
  REF_INT *read_id;
  REF_DBL *read_param;
  REF_DBL filler;

  REF_INT ngeom_read, ngeom_keep;
  REF_INT section_size;

  REF_INT *dest;
  REF_INT *geom_to_send;
  REF_INT *start_to_send;

  REF_INT geom_to_receive;

  REF_INT geom, i;
  REF_INT part, node;
  REF_INT new_location;

  chunk = MAX(1000000, ngeom / ref_mpi_n(ref_mpi));
  chunk = MIN(chunk, ngeom);

  ref_malloc(sent_node, chunk, REF_INT);
  ref_malloc(sent_id, chunk, REF_INT);
  ref_malloc(sent_param, 2 * chunk, REF_DBL);

  if (ref_mpi_once(ref_mpi)) {
    ref_malloc(geom_to_send, ref_mpi_n(ref_mpi), REF_INT);
    ref_malloc(start_to_send, ref_mpi_n(ref_mpi), REF_INT);
    ref_malloc(read_node, chunk, REF_INT);
    ref_malloc(read_id, chunk, REF_INT);
    ref_malloc(read_param, 2 * chunk, REF_DBL);
    ref_malloc(dest, chunk, REF_INT);

    ngeom_read = 0;
    while (ngeom_read < ngeom) {
      section_size = MIN(chunk, ngeom - ngeom_read);
      for (geom = 0; geom < section_size; geom++) {
        REIS(1, fread(&(read_node[geom]), sizeof(REF_INT), 1, file), "n");
        REIS(1, fread(&(read_id[geom]), sizeof(REF_INT), 1, file), "n");
        for (i = 0; i < 2; i++)
          read_param[i + 2 * geom] = 0.0; /* ensure init */
        for (i = 0; i < type; i++)
          REIS(1, fread(&(read_param[i + 2 * geom]), sizeof(double), 1, file),
               "param");
        if (0 < type)
          REIS(1, fread(&(filler), sizeof(double), 1, file), "filler");
      }
      for (geom = 0; geom < section_size; geom++) read_node[geom]--;

      ngeom_read += section_size;

      for (geom = 0; geom < section_size; geom++)
        dest[geom] =
            ref_part_implicit(nnode, ref_mpi_n(ref_mpi), read_node[geom]);

      each_ref_mpi_part(ref_mpi, part) geom_to_send[part] = 0;
      for (geom = 0; geom < section_size; geom++) geom_to_send[dest[geom]]++;

      start_to_send[0] = 0;
      each_ref_mpi_worker(ref_mpi, part) start_to_send[part] =
          start_to_send[part - 1] + geom_to_send[part - 1];

      each_ref_mpi_part(ref_mpi, part) geom_to_send[part] = 0;
      for (geom = 0; geom < section_size; geom++) {
        new_location = start_to_send[dest[geom]] + geom_to_send[dest[geom]];
        sent_node[new_location] = read_node[geom];
        sent_id[new_location] = read_id[geom];
        sent_param[0 + 2 * new_location] = read_param[0 + 2 * geom];
        sent_param[1 + 2 * new_location] = read_param[1 + 2 * geom];
        geom_to_send[dest[geom]]++;
      }

      /* master keepers */
      ngeom_keep = geom_to_send[0];
      for (geom = 0; geom < ngeom_keep; geom++) {
        RSS(ref_node_local(ref_node, sent_node[geom], &node), "g2l");
        RSS(ref_geom_add(ref_geom, node, type, sent_id[geom],
                         &(sent_param[2 * geom])),
            "add geom");
      }

      /* ship it! */
      each_ref_mpi_worker(ref_mpi, part) if (0 < geom_to_send[part]) {
        RSS(ref_mpi_send(ref_mpi, &(geom_to_send[part]), 1, REF_INT_TYPE, part),
            "send");
        RSS(ref_mpi_send(ref_mpi, &(sent_node[start_to_send[part]]),
                         geom_to_send[part], REF_INT_TYPE, part),
            "send");
        RSS(ref_mpi_send(ref_mpi, &(sent_id[start_to_send[part]]),
                         geom_to_send[part], REF_INT_TYPE, part),
            "send");
        RSS(ref_mpi_send(ref_mpi, &(sent_param[2 * start_to_send[part]]),
                         2 * geom_to_send[part], REF_DBL_TYPE, part),
            "send");
      }
    }

    ref_free(dest);
    ref_free(read_param);
    ref_free(read_id);
    ref_free(read_node);
    ref_free(start_to_send);
    ref_free(geom_to_send);

    /* signal we are done */
    each_ref_mpi_worker(ref_mpi, part) RSS(
        ref_mpi_send(ref_mpi, &end_of_message, 1, REF_INT_TYPE, part), "send");

  } else {
    do {
      RSS(ref_mpi_recv(ref_mpi, &geom_to_receive, 1, REF_INT_TYPE, 0), "recv");
      if (geom_to_receive > 0) {
        RSS(ref_mpi_recv(ref_mpi, sent_node, geom_to_receive, REF_INT_TYPE, 0),
            "send");
        RSS(ref_mpi_recv(ref_mpi, sent_id, geom_to_receive, REF_INT_TYPE, 0),
            "send");
        RSS(ref_mpi_recv(ref_mpi, sent_param, 2 * geom_to_receive, REF_DBL_TYPE,
                         0),
            "send");
        for (geom = 0; geom < geom_to_receive; geom++) {
          RSS(ref_node_local(ref_node, sent_node[geom], &node), "g2l");
          RSS(ref_geom_add(ref_geom, node, type, sent_id[geom],
                           &(sent_param[2 * geom])),
              "add geom");
        }
      }
    } while (geom_to_receive != end_of_message);
  }

  free(sent_param);
  free(sent_id);
  free(sent_node);

  return REF_SUCCESS;
}

static REF_STATUS ref_part_meshb_geom_bcast(REF_GEOM ref_geom, REF_INT ngeom,
                                            REF_INT type, REF_NODE ref_node,
                                            FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT chunk;
  REF_INT *read_node;
  REF_INT *read_id;
  REF_DBL *read_param;
  REF_DBL filler;

  REF_INT ngeom_read;
  REF_INT section_size;

  REF_INT geom;
  REF_INT i, local;

  chunk = MAX(1000000, ngeom / ref_mpi_n(ref_mpi));
  chunk = MIN(chunk, ngeom);

  ref_malloc(read_node, chunk, REF_INT);
  ref_malloc(read_id, chunk, REF_INT);
  ref_malloc(read_param, 2 * chunk, REF_DBL);

  ngeom_read = 0;
  while (ngeom_read < ngeom) {
    section_size = MIN(chunk, ngeom - ngeom_read);
    if (ref_mpi_once(ref_mpi)) {
      for (geom = 0; geom < section_size; geom++) {
        REIS(1, fread(&(read_node[geom]), sizeof(REF_INT), 1, file), "n");
        REIS(1, fread(&(read_id[geom]), sizeof(REF_INT), 1, file), "n");
        for (i = 0; i < 2; i++)
          read_param[i + 2 * geom] = 0.0; /* ensure init */
        for (i = 0; i < type; i++)
          REIS(1, fread(&(read_param[i + 2 * geom]), sizeof(double), 1, file),
               "param");
        if (0 < type)
          REIS(1, fread(&(filler), sizeof(double), 1, file), "filler");
      }
      for (geom = 0; geom < section_size; geom++) read_node[geom]--;
    }
    RSS(ref_mpi_bcast(ref_mpi, read_node, section_size, REF_INT_TYPE), "nd");
    RSS(ref_mpi_bcast(ref_mpi, read_id, section_size, REF_INT_TYPE), "id");
    RSS(ref_mpi_bcast(ref_mpi, read_param, 2 * section_size, REF_DBL_TYPE),
        "pm");
    for (geom = 0; geom < section_size; geom++) {
      RXS(ref_node_local(ref_node, read_node[geom], &local), REF_NOT_FOUND,
          "local");
      if (REF_EMPTY != local) {
        RSS(ref_geom_add(ref_geom, local, type, read_id[geom],
                         &(read_param[2 * geom])),
            "add geom");
      }
    }
    ngeom_read += section_size;
  }

  ref_free(read_param);
  ref_free(read_id);
  ref_free(read_node);

  return REF_SUCCESS;
}

static REF_STATUS ref_part_meshb_cell(REF_CELL ref_cell, REF_INT ncell,
                                      REF_NODE ref_node, REF_INT nnode,
                                      FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT ncell_read;
  REF_INT chunk;
  REF_INT end_of_message = REF_EMPTY;
  REF_INT elements_to_receive;
  REF_INT *c2n;
  REF_INT *c2t;
  REF_GLOB *sent_c2n;
  REF_INT *dest;
  REF_INT *sent_part;
  REF_INT *elements_to_send;
  REF_INT *start_to_send;
  REF_INT node_per, size_per;
  REF_INT section_size;
  REF_INT cell;
  REF_INT part, node;
  REF_INT ncell_keep;
  REF_INT new_location;

  chunk = MAX(1000000, ncell / ref_mpi_n(ref_mpi));

  size_per = ref_cell_size_per(ref_cell);
  node_per = ref_cell_node_per(ref_cell);

  ref_malloc(sent_c2n, size_per * chunk, REF_GLOB);

  if (ref_mpi_once(ref_mpi)) {
    ref_malloc(elements_to_send, ref_mpi_n(ref_mpi), REF_INT);
    ref_malloc(start_to_send, ref_mpi_n(ref_mpi), REF_INT);
    ref_malloc(c2n, size_per * chunk, REF_INT);
    ref_malloc(c2t, (node_per + 1) * chunk, REF_INT);
    ref_malloc(dest, chunk, REF_INT);

    ncell_read = 0;
    while (ncell_read < ncell) {
      section_size = MIN(chunk, ncell - ncell_read);
      if (node_per == size_per) {
        RES((size_t)(section_size * (node_per + 1)),
            fread(c2t, sizeof(REF_INT), (size_t)(section_size * (node_per + 1)),
                  file),
            "cn");
        for (cell = 0; cell < section_size; cell++)
          for (node = 0; node < node_per; node++)
            c2n[node + size_per * cell] = c2t[node + (node_per + 1) * cell];
      } else {
        RES((size_t)(section_size * size_per),
            fread(c2n, sizeof(REF_INT), (size_t)(section_size * size_per),
                  file),
            "cn");
      }
      for (cell = 0; cell < section_size; cell++)
        for (node = 0; node < node_per; node++) c2n[node + size_per * cell]--;

      ncell_read += section_size;

      for (cell = 0; cell < section_size; cell++)
        dest[cell] =
            ref_part_implicit(nnode, ref_mpi_n(ref_mpi), c2n[size_per * cell]);

      each_ref_mpi_part(ref_mpi, part) elements_to_send[part] = 0;
      for (cell = 0; cell < section_size; cell++)
        elements_to_send[dest[cell]]++;

      start_to_send[0] = 0;
      each_ref_mpi_worker(ref_mpi, part) start_to_send[part] =
          start_to_send[part - 1] + elements_to_send[part - 1];

      each_ref_mpi_part(ref_mpi, part) elements_to_send[part] = 0;
      for (cell = 0; cell < section_size; cell++) {
        new_location = start_to_send[dest[cell]] + elements_to_send[dest[cell]];
        for (node = 0; node < size_per; node++)
          sent_c2n[node + size_per * new_location] =
              c2n[node + size_per * cell];
        elements_to_send[dest[cell]]++;
      }

      /* master keepers */

      ncell_keep = elements_to_send[0];
      if (0 < ncell_keep) {
        ref_malloc_init(sent_part, size_per * ncell_keep, REF_INT, REF_EMPTY);

        for (cell = 0; cell < ncell_keep; cell++)
          for (node = 0; node < node_per; node++)
            sent_part[node + size_per * cell] = ref_part_implicit(
                nnode, ref_mpi_n(ref_mpi), sent_c2n[node + size_per * cell]);

        RSS(ref_cell_add_many_global(ref_cell, ref_node, ncell_keep, sent_c2n,
                                     sent_part, ref_mpi_rank(ref_mpi)),
            "glob");

        ref_free(sent_part);
      }

      each_ref_mpi_worker(ref_mpi, part) if (0 < elements_to_send[part]) {
        RSS(ref_mpi_send(ref_mpi, &(elements_to_send[part]), 1, REF_INT_TYPE,
                         part),
            "send");
        RSS(ref_mpi_send(ref_mpi, &(sent_c2n[size_per * start_to_send[part]]),
                         size_per * elements_to_send[part], REF_GLOB_TYPE,
                         part),
            "send");
      }
    }

    ref_free(dest);
    ref_free(c2t);
    ref_free(c2n);
    ref_free(start_to_send);
    ref_free(elements_to_send);

    /* signal we are done */
    each_ref_mpi_worker(ref_mpi, part) RSS(
        ref_mpi_send(ref_mpi, &end_of_message, 1, REF_INT_TYPE, part), "send");

  } else {
    do {
      RSS(ref_mpi_recv(ref_mpi, &elements_to_receive, 1, REF_INT_TYPE, 0),
          "recv");
      if (elements_to_receive > 0) {
        RSS(ref_mpi_recv(ref_mpi, sent_c2n, size_per * elements_to_receive,
                         REF_GLOB_TYPE, 0),
            "send");

        ref_malloc_init(sent_part, size_per * elements_to_receive, REF_INT,
                        REF_EMPTY);

        for (cell = 0; cell < elements_to_receive; cell++)
          for (node = 0; node < node_per; node++)
            sent_part[node + size_per * cell] = ref_part_implicit(
                nnode, ref_mpi_n(ref_mpi), sent_c2n[node + size_per * cell]);

        RSS(ref_cell_add_many_global(ref_cell, ref_node, elements_to_receive,
                                     sent_c2n, sent_part,
                                     ref_mpi_rank(ref_mpi)),
            "many glob");

        ref_free(sent_part);
      }
    } while (elements_to_receive != end_of_message);
  }

  free(sent_c2n);

  RSS(ref_migrate_shufflin_cell(ref_node, ref_cell), "fill ghosts");

  return REF_SUCCESS;
}

static REF_STATUS ref_part_meshb_cell_bcast(REF_CELL ref_cell, REF_INT ncell,
                                            REF_NODE ref_node, FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT ncell_read;
  REF_INT chunk;
  REF_INT section_size;
  REF_INT *c2n;
  REF_INT *c2t;
  REF_INT node_per, size_per;
  REF_INT cell, node, local, new_cell;
  REF_BOOL have_all_nodes, one_node_local;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  chunk = MAX(1000000, ncell / ref_mpi_n(ref_mpi));

  size_per = ref_cell_size_per(ref_cell);
  node_per = ref_cell_node_per(ref_cell);

  ref_malloc(c2n, size_per * chunk, REF_INT);

  ncell_read = 0;
  while (ncell_read < ncell) {
    section_size = MIN(chunk, ncell - ncell_read);
    if (ref_mpi_once(ref_mpi)) {
      if (node_per == size_per) {
        ref_malloc(c2t, (node_per + 1) * section_size, REF_INT);
        RES((size_t)(section_size * (node_per + 1)),
            fread(c2t, sizeof(REF_INT), (size_t)(section_size * (node_per + 1)),
                  file),
            "cn");
        for (cell = 0; cell < section_size; cell++)
          for (node = 0; node < node_per; node++)
            c2n[node + size_per * cell] = c2t[node + (node_per + 1) * cell];
        ref_free(c2t);
      } else {
        RES((size_t)(section_size * size_per),
            fread(c2n, sizeof(REF_INT), (size_t)(section_size * size_per),
                  file),
            "cn");
      }
      for (cell = 0; cell < section_size; cell++)
        for (node = 0; node < node_per; node++) c2n[node + size_per * cell]--;
    }
    RSS(ref_mpi_bcast(ref_mpi, c2n, size_per * section_size, REF_INT_TYPE),
        "broadcast read c2n");

    /* convert to local nodes and add if local */
    for (cell = 0; cell < section_size; cell++) {
      have_all_nodes = REF_TRUE;
      for (node = 0; node < node_per; node++) {
        RXS(ref_node_local(ref_node, c2n[node + size_per * cell], &local),
            REF_NOT_FOUND, "local");
        if (REF_EMPTY != local) {
          nodes[node] = local;
        } else {
          have_all_nodes = REF_FALSE;
          break;
        }
      }
      if (have_all_nodes) {
        if (node_per != size_per)
          nodes[node_per] = c2n[node_per + size_per * cell];
        one_node_local = REF_FALSE;
        for (node = 0; node < node_per; node++) {
          one_node_local =
              (one_node_local || ref_node_owned(ref_node, nodes[node]));
        }
        if (one_node_local)
          RSS(ref_cell_add(ref_cell, nodes, &new_cell), "add");
      }
    }

    ncell_read += section_size;
  }

  free(c2n);

  return REF_SUCCESS;
}

static REF_STATUS ref_part_meshb(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                 const char *filename) {
  REF_BOOL verbose = REF_FALSE;
  REF_INT version, dim;
  REF_BOOL available;
  REF_FILEPOS next_position;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_GEOM ref_geom;
  FILE *file;
  REF_BOOL swap_endian = REF_FALSE;
  REF_BOOL has_id = REF_TRUE;
  REF_INT nnode, ncell;
  REF_INT type, geom_keyword, ngeom;
  REF_INT cad_data_keyword;

  file = NULL;
  if (ref_mpi_once(ref_mpi)) {
    RSS(ref_import_meshb_header(filename, &version, key_pos), "header");
    if (verbose) printf("meshb version %d\n", version);
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
  }
  RSS(ref_mpi_bcast(ref_mpi, &dim, 1, REF_INT_TYPE), "bcast");
  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create grid");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);
  ref_geom = ref_grid_geom(ref_grid);
  ref_grid_twod(ref_grid) = (2 == dim);

  if (ref_grid_once(ref_grid)) {
    RSS(ref_import_meshb_jump(file, version, key_pos, 4, &available,
                              &next_position),
        "jump");
    RAS(available, "meshb missing vertex");
    REIS(1, fread((unsigned char *)&nnode, 4, 1, file), "nnode");
    if (verbose) printf("nnode %d\n", nnode);
  }
  RSS(ref_mpi_bcast(ref_mpi, &nnode, 1, REF_INT_TYPE), "bcast");
  RSS(ref_part_node(file, swap_endian, has_id, ref_grid_twod(ref_grid),
                    ref_node, nnode),
      "part node");
  if (ref_grid_once(ref_grid))
    REIS(next_position, ftello(file), "end location");

  if (ref_grid_once(ref_grid)) {
    RSS(ref_import_meshb_jump(file, version, key_pos, 8, &available,
                              &next_position),
        "jump");
    if (available) {
      REIS(1, fread((unsigned char *)&ncell, 4, 1, file), "ntet");
      if (verbose) printf("ntet %d\n", ncell);
    }
  }
  RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");
  if (available) {
    RSS(ref_mpi_bcast(ref_mpi, &ncell, 1, REF_INT_TYPE), "bcast");
    RSS(ref_part_meshb_cell(ref_grid_tet(ref_grid), ncell, ref_node, nnode,
                            file),
        "part cell");
    if (ref_grid_once(ref_grid))
      REIS(next_position, ftello(file), "end location");
  }

  if (ref_grid_once(ref_grid)) {
    RSS(ref_import_meshb_jump(file, version, key_pos, 6, &available,
                              &next_position),
        "jump");
    if (available) {
      REIS(1, fread((unsigned char *)&ncell, 4, 1, file), "ntri");
      if (verbose) printf("ntri %d\n", ncell);
    }
  }
  RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");
  if (available) {
    RSS(ref_mpi_bcast(ref_mpi, &ncell, 1, REF_INT_TYPE), "bcast");
    RSS(ref_part_meshb_cell(ref_grid_tri(ref_grid), ncell, ref_node, nnode,
                            file),
        "part cell");
    if (ref_grid_once(ref_grid))
      REIS(next_position, ftello(file), "end location");
  }

  if (ref_grid_once(ref_grid)) {
    RSS(ref_import_meshb_jump(file, version, key_pos, 5, &available,
                              &next_position),
        "jump");
    if (available) {
      REIS(1, fread((unsigned char *)&ncell, 4, 1, file), "nedge");
      if (verbose) printf("nedge %d\n", ncell);
    }
  }
  RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");
  if (available) {
    RSS(ref_mpi_bcast(ref_mpi, &ncell, 1, REF_INT_TYPE), "bcast");
    RSS(ref_part_meshb_cell(ref_grid_edg(ref_grid), ncell, ref_node, nnode,
                            file),
        "part cell");
    if (ref_grid_once(ref_grid))
      REIS(next_position, ftello(file), "end location");
  }

  each_ref_type(ref_geom, type) {
    if (ref_grid_once(ref_grid)) {
      geom_keyword = 40 + type;
      RSS(ref_import_meshb_jump(file, version, key_pos, geom_keyword,
                                &available, &next_position),
          "jump");
      if (available) {
        REIS(1, fread((unsigned char *)&ngeom, 4, 1, file), "ngeom");
        if (verbose) printf("type %d ngeom %d\n", type, ngeom);
      }
    }
    RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");
    if (available) {
      RSS(ref_mpi_bcast(ref_mpi, &ngeom, 1, REF_INT_TYPE), "bcast");
      RSS(ref_part_meshb_geom(ref_geom, ngeom, type, ref_node, nnode, file),
          "part geom");
      if (ref_grid_once(ref_grid))
        REIS(next_position, ftello(file), "end location");
    }
  }

  if (ref_grid_once(ref_grid)) {
    cad_data_keyword = 126; /* GmfByteFlow */
    RSS(ref_import_meshb_jump(file, version, key_pos, cad_data_keyword,
                              &available, &next_position),
        "jump");
    if (available) {
      REIS(
          1,
          fread((unsigned char *)&ref_geom_cad_data_size(ref_geom), 4, 1, file),
          "cad_data_size");
      if (verbose)
        printf("cad_data_size %ld\n", (long)ref_geom_cad_data_size(ref_geom));
      /* safe non-NULL free, if already allocated, to prevent mem leaks */
      ref_free(ref_geom_cad_data(ref_geom));
      ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                        ref_geom_cad_data_size(ref_geom), REF_BYTE);
      REIS(ref_geom_cad_data_size(ref_geom),
           fread(ref_geom_cad_data(ref_geom), sizeof(REF_BYTE),
                 ref_geom_cad_data_size(ref_geom), file),
           "cad_data");
      REIS(next_position, ftello(file), "end location");
    }
  }
  RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");
  if (available) {
    RSS(ref_mpi_bcast(ref_mpi, &ref_geom_cad_data_size(ref_geom), 1,
                      REF_LONG_TYPE),
        "bcast");
    if (!ref_grid_once(ref_grid))
      ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                        ref_geom_cad_data_size(ref_geom), REF_BYTE);
    RSS(ref_mpi_bcast(ref_mpi, ref_geom_cad_data(ref_geom),
                      (REF_INT)ref_geom_cad_data_size(ref_geom), REF_BYTE_TYPE),
        "bcast");
  }

  RSS(ref_geom_ghost(ref_geom, ref_node), "fill geom ghosts");
  RSS(ref_node_ghost_real(ref_node), "ghost real");

  RSS(ref_grid_inward_boundary_orientation(ref_grid),
      "inward boundary orientation");

  if (ref_grid_once(ref_grid)) {
    fclose(file);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_part_cad_data(REF_GRID ref_grid, const char *filename) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  FILE *file;
  REF_INT version, dim;
  REF_BOOL available;
  REF_FILEPOS next_position;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  REF_DICT ref_dict;
  REF_INT cad_data_keyword;
  REF_BOOL verbose = REF_FALSE;
  size_t end_of_string;

  end_of_string = strlen(filename);

  if (strcmp(&filename[end_of_string - 6], ".meshb") != 0)
    RSS(REF_INVALID, "expected .meshb extension");

  file = NULL;
  if (ref_mpi_once(ref_mpi)) {
    RSS(ref_dict_create(&ref_dict), "create dict");
    RSS(ref_import_meshb_header(filename, &version, key_pos), "header");
    if (verbose) printf("meshb version %d\n", version);
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
  }

  if (ref_grid_once(ref_grid)) {
    cad_data_keyword = 126; /* GmfByteFlow */
    RSS(ref_import_meshb_jump(file, version, key_pos, cad_data_keyword,
                              &available, &next_position),
        "jump");
    if (available) {
      REIS(
          1,
          fread((unsigned char *)&ref_geom_cad_data_size(ref_geom), 4, 1, file),
          "cad_data_size");
      if (verbose)
        printf("cad_data_size %ld\n", (long)ref_geom_cad_data_size(ref_geom));
      /* safe non-NULL free, if already allocated, to prevent mem leaks */
      ref_free(ref_geom_cad_data(ref_geom));
      ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                        ref_geom_cad_data_size(ref_geom), REF_BYTE);
      REIS(ref_geom_cad_data_size(ref_geom),
           fread(ref_geom_cad_data(ref_geom), sizeof(REF_BYTE),
                 (size_t)ref_geom_cad_data_size(ref_geom), file),
           "cad_data");
      REIS(next_position, ftello(file), "end location");
    }
  }
  RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");

  RAS(available, "GmfByteFlow keyword for cad data missing");

  if (available) {
    RSS(ref_mpi_bcast(ref_mpi, &ref_geom_cad_data_size(ref_geom), 1,
                      REF_INT_TYPE),
        "bcast");
    if (!ref_grid_once(ref_grid))
      ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                        ref_geom_cad_data_size(ref_geom), REF_BYTE);
    RSS(ref_mpi_bcast(ref_mpi, ref_geom_cad_data(ref_geom),
                      (REF_INT)ref_geom_cad_data_size(ref_geom), REF_BYTE_TYPE),
        "bcast");
  }

  if (ref_grid_once(ref_grid)) {
    fclose(file);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_part_cad_association(REF_GRID ref_grid, const char *filename) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  FILE *file;
  REF_INT version, dim;
  REF_BOOL available;
  REF_FILEPOS next_position;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  REF_INT type, geom_keyword;
  REF_INT ngeom;
  REF_BOOL verbose = REF_FALSE;
  size_t end_of_string;

  end_of_string = strlen(filename);

  if (strcmp(&filename[end_of_string - 6], ".meshb") != 0)
    RSS(REF_INVALID, "expected .meshb extension");

  file = NULL;
  if (ref_mpi_once(ref_mpi)) {
    RSS(ref_import_meshb_header(filename, &version, key_pos), "header");
    if (verbose) printf("meshb version %d\n", version);
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
  }

  RSS(ref_geom_initialize(ref_geom), "clear out previous assoc");

  each_ref_type(ref_geom, type) {
    if (ref_grid_once(ref_grid)) {
      geom_keyword = 40 + type;
      RSS(ref_import_meshb_jump(file, version, key_pos, geom_keyword,
                                &available, &next_position),
          "jump");
      if (available) {
        REIS(1, fread((unsigned char *)&ngeom, 4, 1, file), "ngeom");
        if (verbose) printf("type %d ngeom %d\n", type, ngeom);
      }
    }
    RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");
    if (available) {
      RSS(ref_mpi_bcast(ref_mpi, &ngeom, 1, REF_INT_TYPE), "bcast");
      RSS(ref_part_meshb_geom_bcast(ref_geom, ngeom, type, ref_node, file),
          "part geom bcast");
      if (ref_grid_once(ref_grid))
        REIS(next_position, ftello(file), "end location");
    }
  }

  if (ref_grid_once(ref_grid)) {
    fclose(file);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_part_cad_discrete_edge(REF_GRID ref_grid, const char *filename) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  FILE *file;
  REF_INT version, dim;
  REF_BOOL available;
  REF_FILEPOS next_position;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  REF_INT ncell;
  REF_BOOL verbose = REF_FALSE;
  size_t end_of_string;

  end_of_string = strlen(filename);

  if (strcmp(&filename[end_of_string - 6], ".meshb") != 0)
    RSS(REF_INVALID, "expected .meshb extension");

  file = NULL;
  if (ref_mpi_once(ref_mpi)) {
    RSS(ref_import_meshb_header(filename, &version, key_pos), "header");
    if (verbose) printf("meshb version %d\n", version);
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
  }

  RSS(ref_cell_free(ref_grid_edg(ref_grid)), "clear out edge");
  RSS(ref_cell_create(&ref_grid_edg(ref_grid), REF_CELL_EDG), "edg");

  if (ref_grid_once(ref_grid)) {
    RSS(ref_import_meshb_jump(file, version, key_pos, 5, &available,
                              &next_position),
        "jump");
    if (available) {
      REIS(1, fread((unsigned char *)&ncell, 4, 1, file), "nedge");
      if (verbose) printf("nedge %d\n", ncell);
    }
  }
  RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");

  RAS(available, "no edge available in meshb");

  RSS(ref_mpi_bcast(ref_mpi, &ncell, 1, REF_INT_TYPE), "bcast");

  RSS(ref_part_meshb_cell_bcast(ref_grid_edg(ref_grid), ncell, ref_node, file),
      "part cell");

  if (ref_grid_once(ref_grid)) {
    REIS(next_position, ftello(file), "end location");
    fclose(file);
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_part_bin_ugrid_pack_cell(
    FILE *file, REF_BOOL swap_endian, REF_BOOL sixty_four_bit,
    REF_FILEPOS conn_offset, REF_FILEPOS faceid_offset, REF_INT section_size,
    REF_LONG ncell_read, REF_INT node_per, REF_INT size_per, REF_GLOB *c2n) {
  REF_INT cell, node;
  REF_FILEPOS ibyte;
  ibyte = (sixty_four_bit ? 8 : 4);

  if (sixty_four_bit) {
    REF_LONG *c2t;
    ref_malloc(c2t, node_per * section_size, REF_LONG);
    REIS(0,
         fseeko(file, conn_offset + ibyte * node_per * (REF_FILEPOS)ncell_read,
                SEEK_SET),
         "seek conn failed");
    RES((size_t)(section_size * node_per),
        fread(c2t, sizeof(REF_LONG), (size_t)(section_size * node_per), file),
        "cn");
    for (cell = 0; cell < section_size * node_per; cell++) {
      if (swap_endian) SWAP_LONG(c2t[cell]);
    }
    for (cell = 0; cell < section_size; cell++) {
      for (node = 0; node < node_per; node++) {
        c2n[node + size_per * cell] = c2t[node + node_per * cell] - 1;
      }
    }
    if (node_per != size_per) {
      REF_LONG *tag;
      ref_malloc(tag, section_size, REF_LONG);
      REIS(0,
           fseeko(file, faceid_offset + ibyte * (REF_FILEPOS)ncell_read,
                  SEEK_SET),
           "seek tag failed");
      RES((size_t)(section_size),
          fread(tag, sizeof(REF_LONG), (size_t)section_size, file), "tag");
      for (cell = 0; cell < section_size; cell++) {
        if (swap_endian) SWAP_LONG(tag[cell]);
      }
      /* sort into right locations */
      for (cell = 0; cell < section_size; cell++)
        c2n[node_per + cell * size_per] = tag[cell];
      ref_free(tag);
    }
    ref_free(c2t);
  } else {
    REF_INT *c2t;
    ref_malloc(c2t, node_per * section_size, REF_INT);
    REIS(0,
         fseeko(file, conn_offset + ibyte * node_per * (REF_FILEPOS)ncell_read,
                SEEK_SET),
         "seek conn failed");
    RES((size_t)(section_size * node_per),
        fread(c2t, sizeof(REF_INT), (size_t)(section_size * node_per), file),
        "cn");
    for (cell = 0; cell < section_size * node_per; cell++) {
      if (swap_endian) SWAP_INT(c2t[cell]);
    }
    for (cell = 0; cell < section_size; cell++) {
      for (node = 0; node < node_per; node++) {
        c2n[node + size_per * cell] = c2t[node + node_per * cell] - 1;
      }
    }
    if (node_per != size_per) {
      REF_INT *tag;
      ref_malloc(tag, section_size, REF_INT);
      REIS(0,
           fseeko(file, faceid_offset + ibyte * (REF_FILEPOS)ncell_read,
                  SEEK_SET),
           "seek tag failed");
      RES((size_t)(section_size),
          fread(tag, sizeof(REF_INT), (size_t)section_size, file), "tag");
      for (cell = 0; cell < section_size; cell++) {
        if (swap_endian) SWAP_INT(tag[cell]);
      }
      /* sort into right locations */
      for (cell = 0; cell < section_size; cell++)
        c2n[node_per + cell * size_per] = tag[cell];
      ref_free(tag);
    }
    ref_free(c2t);
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_part_bin_ugrid_cell(REF_CELL ref_cell, REF_LONG ncell,
                                          REF_NODE ref_node, REF_GLOB nnode,
                                          FILE *file, REF_FILEPOS conn_offset,
                                          REF_FILEPOS faceid_offset,
                                          REF_BOOL swap_endian,
                                          REF_BOOL sixty_four_bit) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_LONG ncell_read;
  REF_INT chunk;
  REF_INT end_of_message = REF_EMPTY;
  REF_INT elements_to_receive;
  REF_GLOB *c2n;
  REF_GLOB *sent_c2n;
  REF_INT *dest;
  REF_INT *sent_part;
  REF_INT *elements_to_send;
  REF_INT *start_to_send;
  REF_INT node_per, size_per;
  REF_INT section_size;
  REF_INT cell;
  REF_INT part, node;
  REF_INT ncell_keep;
  REF_INT new_location;

  chunk = MAX(1000000, (REF_INT)(ncell / ref_mpi_n(ref_node_mpi(ref_node))));

  size_per = ref_cell_size_per(ref_cell);
  node_per = ref_cell_node_per(ref_cell);

  ref_malloc(sent_c2n, size_per * chunk, REF_GLOB);

  if (ref_mpi_once(ref_mpi)) {
    ref_malloc(elements_to_send, ref_mpi_n(ref_mpi), REF_INT);
    ref_malloc(start_to_send, ref_mpi_n(ref_mpi), REF_INT);
    ref_malloc(c2n, size_per * chunk, REF_GLOB);
    ref_malloc(dest, chunk, REF_INT);

    ncell_read = 0;
    while (ncell_read < ncell) {
      section_size = (REF_INT)MIN((REF_LONG)chunk, ncell - ncell_read);

      RSS(ref_part_bin_ugrid_pack_cell(file, swap_endian, sixty_four_bit,
                                       conn_offset, faceid_offset, section_size,
                                       ncell_read, node_per, size_per, c2n),
          "read c2n");

      ncell_read += section_size;

      for (cell = 0; cell < section_size; cell++)
        dest[cell] =
            ref_part_implicit(nnode, ref_mpi_n(ref_mpi), c2n[size_per * cell]);

      each_ref_mpi_part(ref_mpi, part) elements_to_send[part] = 0;
      for (cell = 0; cell < section_size; cell++)
        elements_to_send[dest[cell]]++;

      start_to_send[0] = 0;
      each_ref_mpi_worker(ref_mpi, part) start_to_send[part] =
          start_to_send[part - 1] + elements_to_send[part - 1];

      each_ref_mpi_part(ref_mpi, part) elements_to_send[part] = 0;
      for (cell = 0; cell < section_size; cell++) {
        new_location = start_to_send[dest[cell]] + elements_to_send[dest[cell]];
        for (node = 0; node < size_per; node++)
          sent_c2n[node + size_per * new_location] =
              c2n[node + size_per * cell];
        elements_to_send[dest[cell]]++;
      }

      /* master keepers */

      ncell_keep = elements_to_send[0];
      if (0 < ncell_keep) {
        ref_malloc_init(sent_part, size_per * ncell_keep, REF_INT, REF_EMPTY);

        for (cell = 0; cell < ncell_keep; cell++)
          for (node = 0; node < node_per; node++)
            sent_part[node + size_per * cell] = ref_part_implicit(
                nnode, ref_mpi_n(ref_mpi), sent_c2n[node + size_per * cell]);

        RSS(ref_cell_add_many_global(ref_cell, ref_node, ncell_keep, sent_c2n,
                                     sent_part, ref_mpi_rank(ref_mpi)),
            "glob");

        ref_free(sent_part);
      }

      each_ref_mpi_worker(ref_mpi, part) if (0 < elements_to_send[part]) {
        RSS(ref_mpi_send(ref_mpi, &(elements_to_send[part]), 1, REF_INT_TYPE,
                         part),
            "send");
        RSS(ref_mpi_send(ref_mpi, &(sent_c2n[size_per * start_to_send[part]]),
                         size_per * elements_to_send[part], REF_GLOB_TYPE,
                         part),
            "send");
      }
    }

    ref_free(dest);
    ref_free(c2n);
    ref_free(start_to_send);
    ref_free(elements_to_send);

    /* signal we are done */
    each_ref_mpi_worker(ref_mpi, part) RSS(
        ref_mpi_send(ref_mpi, &end_of_message, 1, REF_INT_TYPE, part), "send");

  } else {
    do {
      RSS(ref_mpi_recv(ref_mpi, &elements_to_receive, 1, REF_INT_TYPE, 0),
          "recv");
      if (elements_to_receive > 0) {
        RSS(ref_mpi_recv(ref_mpi, sent_c2n, size_per * elements_to_receive,
                         REF_GLOB_TYPE, 0),
            "send");

        ref_malloc_init(sent_part, size_per * elements_to_receive, REF_INT,
                        REF_EMPTY);

        for (cell = 0; cell < elements_to_receive; cell++)
          for (node = 0; node < node_per; node++)
            sent_part[node + size_per * cell] = ref_part_implicit(
                nnode, ref_mpi_n(ref_mpi), sent_c2n[node + size_per * cell]);

        RSS(ref_cell_add_many_global(ref_cell, ref_node, elements_to_receive,
                                     sent_c2n, sent_part,
                                     ref_mpi_rank(ref_mpi)),
            "glob");

        ref_free(sent_part);
      }
    } while (elements_to_receive != end_of_message);
  }

  free(sent_c2n);

  RSS(ref_migrate_shufflin_cell(ref_node, ref_cell), "fill ghosts");

  return REF_SUCCESS;
}

static REF_STATUS ref_part_bin_ugrid(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                     const char *filename, REF_BOOL swap_endian,
                                     REF_BOOL sixty_four_bit) {
  FILE *file;
  REF_LONG nnode, ntri, nqua, ntet, npyr, npri, nhex;

  REF_FILEPOS conn_offset, faceid_offset;

  REF_GRID ref_grid;
  REF_NODE ref_node;

  REF_BOOL has_id = REF_FALSE;
  REF_BOOL instrument = REF_FALSE;

  REF_INT single;
  REF_FILEPOS ibyte;
  ibyte = (sixty_four_bit ? 8 : 4);

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create grid");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  if (instrument) ref_mpi_stopwatch_start(ref_grid_mpi(ref_grid));

  /* header */

  file = NULL;
  if (ref_grid_once(ref_grid)) {
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    if (sixty_four_bit) {
      RES(1, fread(&nnode, sizeof(REF_LONG), 1, file), "nnode");
      RES(1, fread(&ntri, sizeof(REF_LONG), 1, file), "ntri");
      RES(1, fread(&nqua, sizeof(REF_LONG), 1, file), "nqua");
      RES(1, fread(&ntet, sizeof(REF_LONG), 1, file), "ntet");
      RES(1, fread(&npyr, sizeof(REF_LONG), 1, file), "npyr");
      RES(1, fread(&npri, sizeof(REF_LONG), 1, file), "npri");
      RES(1, fread(&nhex, sizeof(REF_LONG), 1, file), "nhex");

      if (swap_endian) {
        SWAP_LONG(nnode);
        SWAP_LONG(ntri);
        SWAP_LONG(nqua);
        SWAP_LONG(ntet);
        SWAP_LONG(npyr);
        SWAP_LONG(npri);
        SWAP_LONG(nhex);
      }
    } else {
      RES(1, fread(&single, sizeof(REF_INT), 1, file), "nnode");
      if (swap_endian) SWAP_INT(single);
      nnode = single;
      RES(1, fread(&single, sizeof(REF_INT), 1, file), "ntri");
      if (swap_endian) SWAP_INT(single);
      ntri = single;
      RES(1, fread(&single, sizeof(REF_INT), 1, file), "nqua");
      if (swap_endian) SWAP_INT(single);
      nqua = single;
      RES(1, fread(&single, sizeof(REF_INT), 1, file), "ntet");
      if (swap_endian) SWAP_INT(single);
      ntet = single;
      RES(1, fread(&single, sizeof(REF_INT), 1, file), "npyr");
      if (swap_endian) SWAP_INT(single);
      npyr = single;
      RES(1, fread(&single, sizeof(REF_INT), 1, file), "npri");
      if (swap_endian) SWAP_INT(single);
      npri = single;
      RES(1, fread(&single, sizeof(REF_INT), 1, file), "nhex");
      if (swap_endian) SWAP_INT(single);
      nhex = single;
    }
  }

  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &nnode, 1, REF_LONG_TYPE), "bcast");
  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &ntri, 1, REF_LONG_TYPE), "bcast");
  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &nqua, 1, REF_LONG_TYPE), "bcast");
  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &ntet, 1, REF_LONG_TYPE), "bcast");
  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &npyr, 1, REF_LONG_TYPE), "bcast");
  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &npri, 1, REF_LONG_TYPE), "bcast");
  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &nhex, 1, REF_LONG_TYPE), "bcast");

  RSS(ref_part_node(file, swap_endian, has_id, REF_FALSE, ref_node, nnode),
      "part node");
  if (instrument) ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "nodes");

  if (0 < ntri) {
    conn_offset = 7 * ibyte + (REF_FILEPOS)nnode * (8 * 3);
    faceid_offset = 7 * ibyte + (REF_FILEPOS)nnode * (8 * 3) +
                    (REF_FILEPOS)ntri * 3 * ibyte +
                    (REF_FILEPOS)nqua * 4 * ibyte;
    RSS(ref_part_bin_ugrid_cell(ref_grid_tri(ref_grid), ntri, ref_node, nnode,
                                file, conn_offset, faceid_offset, swap_endian,
                                sixty_four_bit),
        "tri");
  }

  if (0 < nqua) {
    conn_offset = 7 * ibyte + (REF_FILEPOS)nnode * (8 * 3) +
                  (REF_FILEPOS)ntri * 3 * ibyte;
    faceid_offset = 7 * ibyte + (REF_FILEPOS)nnode * (8 * 3) +
                    (REF_FILEPOS)ntri * 4 * ibyte +
                    (REF_FILEPOS)nqua * 4 * ibyte;
    RSS(ref_part_bin_ugrid_cell(ref_grid_qua(ref_grid), nqua, ref_node, nnode,
                                file, conn_offset, faceid_offset, swap_endian,
                                sixty_four_bit),
        "qua");
  }

  if (instrument) ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "bound");

  if (0 < ntet) {
    conn_offset = 7 * ibyte + (REF_FILEPOS)nnode * (8 * 3) +
                  (REF_FILEPOS)ntri * 4 * ibyte + (REF_FILEPOS)nqua * 5 * ibyte;
    faceid_offset = (REF_FILEPOS)REF_EMPTY;
    RSS(ref_part_bin_ugrid_cell(ref_grid_tet(ref_grid), ntet, ref_node, nnode,
                                file, conn_offset, faceid_offset, swap_endian,
                                sixty_four_bit),
        "tet");
  }

  if (0 < npyr) {
    conn_offset = 7 * ibyte + (REF_FILEPOS)nnode * (8 * 3) +
                  (REF_FILEPOS)ntri * 4 * ibyte +
                  (REF_FILEPOS)nqua * 5 * ibyte + (REF_FILEPOS)ntet * 4 * ibyte;
    faceid_offset = (REF_FILEPOS)REF_EMPTY;
    RSS(ref_part_bin_ugrid_cell(ref_grid_pyr(ref_grid), npyr, ref_node, nnode,
                                file, conn_offset, faceid_offset, swap_endian,
                                sixty_four_bit),
        "pyr");
  }

  if (0 < npri) {
    conn_offset = 7 * ibyte + (REF_FILEPOS)nnode * (8 * 3) +
                  (REF_FILEPOS)ntri * 4 * ibyte +
                  (REF_FILEPOS)nqua * 5 * ibyte +
                  (REF_FILEPOS)ntet * 4 * ibyte + (REF_FILEPOS)npyr * 5 * ibyte;
    faceid_offset = (REF_FILEPOS)REF_EMPTY;
    RSS(ref_part_bin_ugrid_cell(ref_grid_pri(ref_grid), npri, ref_node, nnode,
                                file, conn_offset, faceid_offset, swap_endian,
                                sixty_four_bit),
        "pri");
  }

  if (0 < nhex) {
    conn_offset = 7 * ibyte + (REF_FILEPOS)nnode * (8 * 3) +
                  (REF_FILEPOS)ntri * 4 * ibyte +
                  (REF_FILEPOS)nqua * 5 * ibyte +
                  (REF_FILEPOS)ntet * 4 * ibyte +
                  (REF_FILEPOS)npyr * 5 * ibyte + (REF_FILEPOS)npri * 6 * ibyte;
    faceid_offset = REF_EMPTY;
    RSS(ref_part_bin_ugrid_cell(ref_grid_hex(ref_grid), nhex, ref_node, nnode,
                                file, conn_offset, faceid_offset, swap_endian,
                                sixty_four_bit),
        "hex");
  }

  if (ref_grid_once(ref_grid)) REIS(0, fclose(file), "close file");

  /* ghost xyz */

  RSS(ref_node_ghost_real(ref_node), "ghost real");

  RSS(ref_grid_inward_boundary_orientation(ref_grid),
      "inward boundary orientation");

  if (instrument) ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "volume");

  return REF_SUCCESS;
}

REF_STATUS ref_part_metric_solb(REF_NODE ref_node, const char *filename) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_FILEPOS next_position;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  FILE *file;
  REF_INT chunk;
  REF_DBL *metric;
  REF_INT nnode_read, section_size;
  REF_INT node, local, global;
  REF_BOOL available;
  REF_INT version, dim, nnode, ntype, type;

  file = NULL;
  if (ref_mpi_once(ref_mpi)) {
    RSS(ref_import_meshb_header(filename, &version, key_pos), "header");
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");
    RSS(ref_import_meshb_jump(file, version, key_pos, 3, &available,
                              &next_position),
        "jump");
    RAS(available, "meshb missing dimension");
    REIS(1, fread((unsigned char *)&dim, 4, 1, file), "dim");
    if (dim < 2 || 3 < dim) {
      printf("dim %d not supported\n", dim);
      THROW("dim");
    }

    RSS(ref_import_meshb_jump(file, version, key_pos, 62, &available,
                              &next_position),
        "jump");
    RAS(available, "SolAtVertices missing");
    REIS(1, fread((unsigned char *)&nnode, 4, 1, file), "nnode");
    REIS(1, fread((unsigned char *)&ntype, 4, 1, file), "ntype");
    REIS(1, fread((unsigned char *)&type, 4, 1, file), "type");
    if ((nnode != ref_node_n_global(ref_node)) &&
        (nnode / 2 != ref_node_n_global(ref_node))) {
      printf("file %d ref_node " REF_GLOB_FMT "\n", nnode,
             ref_node_n_global(ref_node));
      THROW("global count mismatch");
    }
    REIS(1, ntype, "number of solutions");
    REIS(3, type, "metric solution type");
  }
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), &dim, 1, REF_INT_TYPE),
      "bcast dim");
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), &nnode, 1, REF_INT_TYPE),
      "bcast nnode");

  chunk = MAX(100000, nnode / ref_mpi_n(ref_node_mpi(ref_node)));
  chunk = MIN(chunk, nnode);

  ref_malloc_init(metric, 6 * chunk, REF_DBL, -1.0);

  nnode_read = 0;
  while (nnode_read < nnode) {
    section_size = MIN(chunk, nnode - nnode_read);
    if (ref_mpi_once(ref_node_mpi(ref_node))) {
      for (node = 0; node < section_size; node++) {
        if (3 == dim) {
          REIS(1, fread(&(metric[0 + 6 * node]), sizeof(REF_DBL), 1, file),
               "m11");
          REIS(1, fread(&(metric[1 + 6 * node]), sizeof(REF_DBL), 1, file),
               "m12");
          /* transposed 3,2 */
          REIS(1, fread(&(metric[3 + 6 * node]), sizeof(REF_DBL), 1, file),
               "m22");
          REIS(1, fread(&(metric[2 + 6 * node]), sizeof(REF_DBL), 1, file),
               "m31");
          REIS(1, fread(&(metric[4 + 6 * node]), sizeof(REF_DBL), 1, file),
               "m32");
          REIS(1, fread(&(metric[5 + 6 * node]), sizeof(REF_DBL), 1, file),
               "m33");
        } else {
          REIS(1, fread(&(metric[0 + 6 * node]), sizeof(REF_DBL), 1, file),
               "m11");
          metric[1 + 6 * node] = 0.0; /* m12 */
          REIS(1, fread(&(metric[2 + 6 * node]), sizeof(REF_DBL), 1, file),
               "m31");
          metric[3 + 6 * node] = 1.0; /* m22 */
          metric[4 + 6 * node] = 0.0; /* m32 */
          REIS(1, fread(&(metric[5 + 6 * node]), sizeof(REF_DBL), 1, file),
               "m33");
        }
      }
      RSS(ref_mpi_bcast(ref_node_mpi(ref_node), metric, 6 * chunk,
                        REF_DBL_TYPE),
          "bcast");
    } else {
      RSS(ref_mpi_bcast(ref_node_mpi(ref_node), metric, 6 * chunk,
                        REF_DBL_TYPE),
          "bcast");
    }
    for (node = 0; node < section_size; node++) {
      global = node + nnode_read;
      RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
      if (REF_EMPTY != local) {
        RSS(ref_node_metric_set(ref_node, local, &(metric[6 * node])),
            "set local node met");
      }
      if (2 == dim) {
        global = nnode + node + nnode_read;
        RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
        if (REF_EMPTY != local) {
          RSS(ref_node_metric_set(ref_node, local, &(metric[6 * node])),
              "set local node met");
        }
      }
    }
    nnode_read += section_size;
  }

  ref_free(metric);
  if (ref_mpi_once(ref_mpi)) {
    REIS(next_position, ftello(file), "end location");
    REIS(0, fclose(file), "close file");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_part_metric(REF_NODE ref_node, const char *filename) {
  FILE *file;
  REF_INT chunk;
  REF_DBL *metric;
  REF_INT section_size;
  REF_GLOB nnode_read, global;
  REF_INT node, local;
  size_t end_of_string;
  REF_BOOL sol_format, found_keyword;
  REF_INT nnode, ntype, type;
  REF_INT status;
  char line[1024];
  REF_BOOL solb_format = REF_FALSE;

  if (ref_mpi_once(ref_node_mpi(ref_node))) {
    end_of_string = strlen(filename);
    if (strcmp(&filename[end_of_string - 5], ".solb") == 0)
      solb_format = REF_TRUE;
  }
  RSS(ref_mpi_all_or(ref_node_mpi(ref_node), &solb_format), "bcast");

  if (solb_format) {
    RSS(ref_part_metric_solb(ref_node, filename), "-metric.solb");
    return REF_SUCCESS;
  }

  file = NULL;
  sol_format = REF_FALSE;
  if (ref_mpi_once(ref_node_mpi(ref_node))) {
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    end_of_string = strlen(filename);
    if (strcmp(&filename[end_of_string - 4], ".sol") == 0) {
      sol_format = REF_TRUE;
      found_keyword = REF_FALSE;
      while (!feof(file)) {
        status = fscanf(file, "%s", line);
        if (EOF == status) break;
        REIS(1, status, "line read failed");

        if (0 == strcmp("SolAtVertices", line)) {
          REIS(1, fscanf(file, "%d", &nnode), "read nnode");
          REIS(ref_node_n_global(ref_node), nnode,
               "wrong vertex number in .sol");
          REIS(2, fscanf(file, "%d %d", &ntype, &type), "read header");
          REIS(1, ntype, "expected one type in .sol");
          REIS(3, type, "expected type GmfSymMat in .sol");
          RAS(0 <= fscanf(file, "%*[^1234567890-+.]"), "skip blank line");
          found_keyword = REF_TRUE;
          break;
        }
      }
      RAS(found_keyword, "SolAtVertices keyword missing from .sol metric");
    }
  }

  chunk = (REF_INT)MAX(
      100000, ref_node_n_global(ref_node) / ref_mpi_n(ref_node_mpi(ref_node)));
  chunk = (REF_INT)MIN((REF_GLOB)chunk, ref_node_n_global(ref_node));

  ref_malloc_init(metric, 6 * chunk, REF_DBL, -1.0);

  nnode_read = 0;
  while (nnode_read < ref_node_n_global(ref_node)) {
    section_size =
        (REF_INT)MIN((REF_GLOB)chunk, ref_node_n_global(ref_node) - nnode_read);
    if (ref_mpi_once(ref_node_mpi(ref_node))) {
      for (node = 0; node < section_size; node++)
        if (sol_format) {
          REIS(6,
               fscanf(file, "%lf %lf %lf %lf %lf %lf", &(metric[0 + 6 * node]),
                      &(metric[1 + 6 * node]),
                      &(metric[3 + 6 * node]), /* transposed 3,2 */
                      &(metric[2 + 6 * node]), &(metric[4 + 6 * node]),
                      &(metric[5 + 6 * node])),
               "metric read error");
        } else {
          REIS(6,
               fscanf(file, "%lf %lf %lf %lf %lf %lf", &(metric[0 + 6 * node]),
                      &(metric[1 + 6 * node]), &(metric[2 + 6 * node]),
                      &(metric[3 + 6 * node]), &(metric[4 + 6 * node]),
                      &(metric[5 + 6 * node])),
               "metric read error");
        }
      RSS(ref_mpi_bcast(ref_node_mpi(ref_node), metric, 6 * chunk,
                        REF_DBL_TYPE),
          "bcast");
    } else {
      RSS(ref_mpi_bcast(ref_node_mpi(ref_node), metric, 6 * chunk,
                        REF_DBL_TYPE),
          "bcast");
    }
    for (node = 0; node < section_size; node++) {
      global = node + nnode_read;
      RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
      if (REF_EMPTY != local) {
        RSS(ref_node_metric_set(ref_node, local, &(metric[6 * node])),
            "set local node met");
      }
    }
    nnode_read += section_size;
  }

  ref_free(metric);
  if (ref_mpi_once(ref_node_mpi(ref_node))) REIS(0, fclose(file), "close file");

  return REF_SUCCESS;
}

REF_STATUS ref_part_bamg_metric(REF_GRID ref_grid, const char *filename) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  FILE *file;
  REF_INT chunk;
  REF_DBL *metric;
  REF_INT file_nnode, section_size;
  REF_GLOB nnode, nnode_read, global;
  REF_INT node, local;
  REF_INT nterm;

  nnode = ref_node_n_global(ref_node) / 2;

  file = NULL;
  if (ref_grid_once(ref_grid)) {
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");
    REIS(2, fscanf(file, "%d %d", &file_nnode, &nterm), "read header");
    REIS(nnode, file_nnode, "wrong node count");
    REIS(3, nterm, "expected 3 term 2x2 anisotropic M");
  }

  chunk = (REF_INT)MAX(100000, nnode / ref_mpi_n(ref_grid_mpi(ref_grid)));
  chunk = (REF_INT)MIN((REF_GLOB)chunk, nnode);

  ref_malloc_init(metric, 3 * chunk, REF_DBL, -1.0);

  nnode_read = 0;
  while (nnode_read < nnode) {
    section_size =
        (REF_INT)MIN((REF_GLOB)chunk, ref_node_n_global(ref_node) - nnode_read);
    if (ref_grid_once(ref_grid)) {
      for (node = 0; node < section_size; node++) {
        REIS(3,
             fscanf(file, "%lf %lf %lf", &(metric[0 + 3 * node]),
                    &(metric[1 + 3 * node]), &(metric[2 + 3 * node])),
             "metric read error");
      }
      RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), metric, 3 * chunk,
                        REF_DBL_TYPE),
          "bcast");
    } else {
      RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), metric, 3 * chunk,
                        REF_DBL_TYPE),
          "bcast");
    }
    for (node = 0; node < section_size; node++) {
      global = node + nnode_read;
      RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
      if (REF_EMPTY != local) {
        RSS(ref_node_metric_form(ref_node, local, metric[0 + 3 * node],
                                 metric[1 + 3 * node], 0, metric[2 + 3 * node],
                                 0, 1),
            "set local node met");
      }
    }
    nnode_read += section_size;
  }

  ref_free(metric);
  if (ref_grid_once(ref_grid)) REIS(0, fclose(file), "close file");

  return REF_SUCCESS;
}

static REF_STATUS ref_part_scalar_solb(REF_NODE ref_node, REF_INT *ldim,
                                       REF_DBL **scalar, const char *filename) {
  REF_FILEPOS next_position;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  FILE *file;
  REF_INT chunk;
  REF_DBL *data;
  REF_INT section_size;
  REF_GLOB nnode_read, global;
  REF_INT node, local;
  size_t end_of_string;
  REF_BOOL available;
  REF_INT version, dim, nnode, ntype, type, i;

  file = NULL;
  if (ref_mpi_once(ref_node_mpi(ref_node))) {
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    end_of_string = strlen(filename);
    REIS(0, strcmp(&filename[end_of_string - 5], ".solb"),
         "solb extension expected");

    RSS(ref_import_meshb_header(filename, &version, key_pos), "head");
    RSS(ref_import_meshb_jump(file, version, key_pos, 3, &available,
                              &next_position),
        "jump");
    RAS(available, "solb missing dimension");
    REIS(1, fread((unsigned char *)&dim, 4, 1, file), "dim");
    REIS(3, dim, "only 3D supported");

    RSS(ref_import_meshb_jump(file, version, key_pos, 62, &available,
                              &next_position),
        "jmp");
    RAS(available, "SolAtVertices missing");
    REIS(1, fread((unsigned char *)&nnode, 4, 1, file), "nnode");
    REIS(1, fread((unsigned char *)&ntype, 4, 1, file), "ntype");
    if ((nnode != ref_node_n_global(ref_node)) &&
        (nnode / 2 != ref_node_n_global(ref_node))) {
      printf("file %d ref_node " REF_GLOB_FMT "\n", nnode,
             ref_node_n_global(ref_node));
      THROW("global count mismatch");
    }
    *ldim = 0;
    for (i = 0; i < ntype; i++) {
      REIS(1, fread((unsigned char *)&type, 4, 1, file), "type");
      if (1 != type && 2 != type) {
        printf("only types 1 (scalar) or 2 (vector) supported for %d type\n",
               type);
      }
      if (1 == type) (*ldim) += 1;
      if (2 == type) (*ldim) += dim;
    }
  }
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), ldim, 1, REF_INT_TYPE),
      "bcast ldim");
  ref_malloc(*scalar, (*ldim) * ref_node_max(ref_node), REF_DBL);

  chunk = (REF_INT)MAX(
      100000, ref_node_n_global(ref_node) / ref_mpi_n(ref_node_mpi(ref_node)));
  chunk = (REF_INT)MIN((REF_GLOB)chunk, ref_node_n_global(ref_node));

  ref_malloc_init(data, (*ldim) * chunk, REF_DBL, -1.0);

  nnode_read = 0;
  while (nnode_read < ref_node_n_global(ref_node)) {
    section_size =
        (REF_INT)MIN((REF_GLOB)chunk, ref_node_n_global(ref_node) - nnode_read);
    if (ref_mpi_once(ref_node_mpi(ref_node))) {
      REIS((*ldim) * section_size,
           fread(data, sizeof(REF_DBL), (size_t)((*ldim) * section_size), file),
           "dat");
      RSS(ref_mpi_bcast(ref_node_mpi(ref_node), data, (*ldim) * chunk,
                        REF_DBL_TYPE),
          "bcast");
    } else {
      RSS(ref_mpi_bcast(ref_node_mpi(ref_node), data, (*ldim) * chunk,
                        REF_DBL_TYPE),
          "bcast");
    }
    for (node = 0; node < section_size; node++) {
      global = node + nnode_read;
      RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
      if (REF_EMPTY != local) {
        for (i = 0; i < *ldim; i++) {
          (*scalar)[i + local * (*ldim)] = data[i + node * (*ldim)];
        }
      }
    }
    nnode_read += section_size;
  }

  ref_free(data);

  if (ref_mpi_once(ref_node_mpi(ref_node))) {
    if (nnode == ref_node_n_global(ref_node))
      REIS(next_position, ftello(file), "end location");
    REIS(0, fclose(file), "close file");
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_part_scalar_snap(REF_NODE ref_node, REF_INT *ldim,
                                       REF_DBL **scalar, const char *filename) {
  REF_FILEPOS next_position;
  FILE *file;
  REF_INT chunk;
  REF_DBL *data;
  REF_INT section_size;
  REF_GLOB nnode_read, global, nnode;
  REF_INT node, local;
  size_t end_of_string;
  REF_INT i;
  unsigned long version_number;
  unsigned long number_of_fields;
  unsigned long number_of_chars;
  char letter;
  unsigned long field_length, uint_nnode;
  int association;

  nnode = -1;
  next_position = -1;
  file = NULL;
  if (ref_mpi_once(ref_node_mpi(ref_node))) {
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    end_of_string = strlen(filename);
    REIS(0, strcmp(&filename[end_of_string - 5], ".snap"),
         "solb extension expected");

    REIS(1, fread(&version_number, sizeof(version_number), 1, file),
         "version number");
    REIS(2, version_number, "only version 2 supported");

    REIS(1, fread(&number_of_fields, sizeof(number_of_fields), 1, file),
         "number");
    REIS(1, number_of_fields, "only one field supported");

    REIS(1, fread(&number_of_chars, sizeof(number_of_chars), 1, file),
         "number");
    for (i = 0; i < (REF_INT)number_of_chars; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "number");
    }

    REIS(1, fread(&field_length, sizeof(field_length), 1, file), "number");
    next_position = (REF_FILEPOS)field_length + ftello(file);
    REIS(1, fread(&uint_nnode, sizeof(uint_nnode), 1, file), "number");
    nnode = (REF_GLOB)uint_nnode;
    REIS(1, fread(&association, sizeof(association), 1, file), "number");

    REIS(-1, association, "field node association only");

    *ldim = 1;

    if ((nnode != ref_node_n_global(ref_node)) &&
        (nnode / 2 != ref_node_n_global(ref_node))) {
      printf("file " REF_GLOB_FMT " ref_node " REF_GLOB_FMT "\n", nnode,
             ref_node_n_global(ref_node));
      THROW("global count mismatch");
    }
  }
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), ldim, 1, REF_INT_TYPE),
      "bcast ldim");
  ref_malloc(*scalar, (*ldim) * ref_node_max(ref_node), REF_DBL);

  chunk = (REF_INT)MAX(
      100000, ref_node_n_global(ref_node) / ref_mpi_n(ref_node_mpi(ref_node)));
  chunk = (REF_INT)MIN((REF_GLOB)chunk, ref_node_n_global(ref_node));

  ref_malloc_init(data, (*ldim) * chunk, REF_DBL, -1.0);

  nnode_read = 0;
  while (nnode_read < ref_node_n_global(ref_node)) {
    section_size =
        (REF_INT)MIN((REF_GLOB)chunk, ref_node_n_global(ref_node) - nnode_read);
    if (ref_mpi_once(ref_node_mpi(ref_node))) {
      REIS((*ldim) * section_size,
           fread(data, sizeof(REF_DBL), (size_t)((*ldim) * section_size), file),
           "dat");
      RSS(ref_mpi_bcast(ref_node_mpi(ref_node), data, (*ldim) * chunk,
                        REF_DBL_TYPE),
          "bcast");
    } else {
      RSS(ref_mpi_bcast(ref_node_mpi(ref_node), data, (*ldim) * chunk,
                        REF_DBL_TYPE),
          "bcast");
    }
    for (node = 0; node < section_size; node++) {
      global = node + nnode_read;
      RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
      if (REF_EMPTY != local) {
        for (i = 0; i < *ldim; i++) {
          (*scalar)[i + local * (*ldim)] = data[i + node * (*ldim)];
        }
      }
    }
    nnode_read += section_size;
  }

  ref_free(data);

  if (ref_mpi_once(ref_node_mpi(ref_node))) {
    if (nnode == ref_node_n_global(ref_node))
      REIS(next_position, ftello(file), "end location");
    REIS(0, fclose(file), "close file");
  }
  return REF_SUCCESS;
}

REF_STATUS ref_part_scalar(REF_NODE ref_node, REF_INT *ldim, REF_DBL **scalar,
                           const char *filename) {
  size_t end_of_string;

  end_of_string = strlen(filename);

  if (end_of_string > 5 && strcmp(&filename[end_of_string - 5], ".solb") == 0) {
    RSS(ref_part_scalar_solb(ref_node, ldim, scalar, filename), "solb failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 5 && strcmp(&filename[end_of_string - 5], ".snap") == 0) {
    RSS(ref_part_scalar_snap(ref_node, ldim, scalar, filename), "snap failed");
    return REF_SUCCESS;
  }

  printf("%s: %d: %s %s\n", __FILE__, __LINE__,
         "input file name extension unknown", filename);
  RSS(REF_FAILURE, "unknown file extension");
  return REF_FAILURE;
}

REF_STATUS ref_part_by_extension(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                 const char *filename) {
  size_t end_of_string;

  end_of_string = strlen(filename);

  if (end_of_string > 10 &&
      strcmp(&filename[end_of_string - 10], ".lb8.ugrid") == 0) {
    RSS(ref_part_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_FALSE,
                           REF_FALSE),
        "lb8_ugrid failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 9 &&
      strcmp(&filename[end_of_string - 9], ".b8.ugrid") == 0) {
    RSS(ref_part_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_TRUE,
                           REF_FALSE),
        "b8_ugrid failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 11 &&
      strcmp(&filename[end_of_string - 11], ".lb8l.ugrid") == 0) {
    RSS(ref_part_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_FALSE,
                           REF_TRUE),
        "lb8_ugrid failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 10 &&
      strcmp(&filename[end_of_string - 10], ".b8l.ugrid") == 0) {
    RSS(ref_part_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_TRUE, REF_TRUE),
        "b8_ugrid failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 12 &&
      strcmp(&filename[end_of_string - 12], ".lb8.ugrid64") == 0) {
    RSS(ref_part_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_FALSE,
                           REF_TRUE),
        "lb8_ugrid failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 11 &&
      strcmp(&filename[end_of_string - 11], ".b8.ugrid64") == 0) {
    RSS(ref_part_bin_ugrid(ref_grid_ptr, ref_mpi, filename, REF_TRUE, REF_TRUE),
        "b8_ugrid failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 6 &&
      strcmp(&filename[end_of_string - 6], ".meshb") == 0) {
    RSS(ref_part_meshb(ref_grid_ptr, ref_mpi, filename), "meshb failed");
    return REF_SUCCESS;
  }
  printf("%s: %d: %s %s\n", __FILE__, __LINE__,
         "input file name extension unknown", filename);
  RSS(REF_FAILURE, "unknown file extension");
  return REF_FAILURE;
}
