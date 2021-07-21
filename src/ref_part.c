
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

#include <math.h>
#include <stdint.h>
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

static REF_STATUS ref_part_meshb_long(FILE *file, REF_INT version,
                                      REF_LONG *value) {
  int int_value;
  long long_value;
  if (version < 4) {
    REIS(1, fread(&int_value, sizeof(int), 1, file), "int value");
    *value = (REF_LONG)int_value;
  } else {
    REIS(1, fread(&long_value, sizeof(long), 1, file), "long value");
    *value = long_value;
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_part_meshb_size(FILE *file, REF_INT version,
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

static REF_STATUS ref_part_node(FILE *file, REF_BOOL swap_endian,
                                REF_INT version, REF_BOOL twod,
                                REF_NODE ref_node, REF_LONG nnode) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT node, new_node;
  REF_INT part;
  REF_INT n;
  REF_LONG id;
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
      if (version > 0) RSS(ref_part_meshb_long(file, version, &id), "nnode");
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
          if (version > 0)
            RSS(ref_part_meshb_long(file, version, &id), "nnode");
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

/* used by SANS? */
REF_STATUS ref_part_meshb_geom_delete_me(REF_GEOM ref_geom, REF_INT ngeom,
                                         REF_INT type, REF_NODE ref_node,
                                         REF_INT nnode, FILE *file);
REF_STATUS ref_part_meshb_geom_delete_me(REF_GEOM ref_geom, REF_INT ngeom,
                                         REF_INT type, REF_NODE ref_node,
                                         REF_INT nnode, FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT end_of_message = REF_EMPTY;
  REF_INT chunk;
  REF_INT *sent_node;
  REF_INT *sent_id;
  REF_INT *sent_gref;
  REF_DBL *sent_param;
  REF_INT *read_node;
  REF_INT *read_id;
  REF_INT *read_gref;
  REF_DBL *read_param;
  double double_gref;

  REF_INT ngeom_read, ngeom_keep;
  REF_INT section_size;

  REF_INT *dest;
  REF_INT *geom_to_send;
  REF_INT *start_to_send;

  REF_INT geom_to_receive;

  REF_INT geom, i, new_geom;
  REF_INT part, node;
  REF_INT new_location;

  chunk = MAX(1000000, ngeom / ref_mpi_n(ref_mpi));
  chunk = MIN(chunk, ngeom);

  ref_malloc(sent_node, chunk, REF_INT);
  ref_malloc(sent_id, chunk, REF_INT);
  ref_malloc(sent_gref, chunk, REF_INT);
  ref_malloc(sent_param, 2 * chunk, REF_DBL);

  if (ref_mpi_once(ref_mpi)) {
    ref_malloc(geom_to_send, ref_mpi_n(ref_mpi), REF_INT);
    ref_malloc(start_to_send, ref_mpi_n(ref_mpi), REF_INT);
    ref_malloc(read_node, chunk, REF_INT);
    ref_malloc(read_id, chunk, REF_INT);
    ref_malloc(read_gref, chunk, REF_INT);
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
        read_gref[geom] = read_id[geom];
        if (0 < type) {
          REIS(1, fread(&(double_gref), sizeof(double), 1, file), "gref");
          read_gref[geom] = (REF_INT)double_gref;
        }
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
        sent_gref[new_location] = read_gref[geom];
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
        RSS(ref_geom_find(ref_geom, node, type, sent_id[geom], &new_geom),
            "find");
        ref_geom_gref(ref_geom, new_geom) = sent_gref[geom];
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
        RSS(ref_mpi_send(ref_mpi, &(sent_gref[start_to_send[part]]),
                         geom_to_send[part], REF_INT_TYPE, part),
            "send");
        RSS(ref_mpi_send(ref_mpi, &(sent_param[2 * start_to_send[part]]),
                         2 * geom_to_send[part], REF_DBL_TYPE, part),
            "send");
      }
    }

    ref_free(dest);
    ref_free(read_param);
    ref_free(read_gref);
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
        RSS(ref_mpi_recv(ref_mpi, sent_gref, geom_to_receive, REF_INT_TYPE, 0),
            "send");
        RSS(ref_mpi_recv(ref_mpi, sent_param, 2 * geom_to_receive, REF_DBL_TYPE,
                         0),
            "send");
        for (geom = 0; geom < geom_to_receive; geom++) {
          RSS(ref_node_local(ref_node, sent_node[geom], &node), "g2l");
          RSS(ref_geom_add(ref_geom, node, type, sent_id[geom],
                           &(sent_param[2 * geom])),
              "add geom");
          RSS(ref_geom_find(ref_geom, node, type, sent_id[geom], &new_geom),
              "find");
          ref_geom_gref(ref_geom, new_geom) = sent_gref[geom];
        }
      }
    } while (geom_to_receive != end_of_message);
  }

  free(sent_param);
  free(sent_gref);
  free(sent_id);
  free(sent_node);

  return REF_SUCCESS;
}

static REF_STATUS ref_part_meshb_geom_bcast(REF_GEOM ref_geom, REF_LONG ngeom,
                                            REF_INT type, REF_NODE ref_node,
                                            REF_INT version, FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_INT chunk;
  REF_LONG *read_node;
  REF_LONG *read_id;
  REF_LONG *read_gref;
  REF_DBL *read_param;
  REF_DBL double_gref;

  REF_LONG ngeom_read;
  REF_INT section_size;

  REF_INT geom, new_geom;
  REF_INT i, local;

  chunk = (REF_INT)MAX(1000000, ngeom / (REF_LONG)ref_mpi_n(ref_mpi));
  chunk = (REF_INT)MIN((REF_LONG)chunk, ngeom);

  ref_malloc(read_node, chunk, REF_LONG);
  ref_malloc(read_id, chunk, REF_LONG);
  ref_malloc(read_gref, chunk, REF_LONG);
  ref_malloc(read_param, 2 * chunk, REF_DBL);

  ngeom_read = 0;
  while (ngeom_read < ngeom) {
    section_size = MIN(chunk, (REF_INT)(ngeom - ngeom_read));
    if (ref_mpi_once(ref_mpi)) {
      for (geom = 0; geom < section_size; geom++) {
        RSS(ref_part_meshb_long(file, version, &(read_node[geom])), "node");
        RSS(ref_part_meshb_long(file, version, &(read_id[geom])), "node");
        for (i = 0; i < 2; i++)
          read_param[i + 2 * geom] = 0.0; /* ensure init */
        for (i = 0; i < type; i++)
          REIS(1, fread(&(read_param[i + 2 * geom]), sizeof(double), 1, file),
               "param");
        read_gref[geom] = read_id[geom];
        if (0 < type) {
          REIS(1, fread(&(double_gref), sizeof(double), 1, file), "gref");
          read_gref[geom] = (REF_LONG)double_gref;
        }
      }
      for (geom = 0; geom < section_size; geom++) read_node[geom]--;
    }
    RSS(ref_mpi_bcast(ref_mpi, read_node, section_size, REF_LONG_TYPE), "nd");
    RSS(ref_mpi_bcast(ref_mpi, read_id, section_size, REF_LONG_TYPE), "id");
    RSS(ref_mpi_bcast(ref_mpi, read_gref, section_size, REF_LONG_TYPE), "gref");
    RSS(ref_mpi_bcast(ref_mpi, read_param, 2 * section_size, REF_DBL_TYPE),
        "pm");
    for (geom = 0; geom < section_size; geom++) {
      RXS(ref_node_local(ref_node, read_node[geom], &local), REF_NOT_FOUND,
          "local");
      if (REF_EMPTY != local) {
        RSS(ref_geom_add(ref_geom, local, type, (REF_INT)read_id[geom],
                         &(read_param[2 * geom])),
            "add geom");
        RSS(ref_geom_find(ref_geom, local, type, (REF_INT)read_id[geom],
                          &new_geom),
            "find");
        ref_geom_gref(ref_geom, new_geom) = (REF_INT)read_gref[geom];
      }
    }
    ngeom_read += section_size;
  }

  ref_free(read_param);
  ref_free(read_gref);
  ref_free(read_id);
  ref_free(read_node);

  return REF_SUCCESS;
}

static REF_STATUS ref_part_meshb_cell(REF_CELL ref_cell, REF_LONG ncell,
                                      REF_NODE ref_node, REF_LONG nnode,
                                      REF_INT version, REF_BOOL pad,
                                      FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_LONG ncell_read;
  REF_INT chunk;
  REF_INT end_of_message = REF_EMPTY;
  REF_INT elements_to_receive;
  REF_GLOB *c2n;
  REF_INT *c2n_int;
  REF_LONG *c2n_long;
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
  size_t nread;

  chunk = (REF_INT)MAX(1000000, ncell / (REF_LONG)ref_mpi_n(ref_mpi));

  node_per = ref_cell_node_per(ref_cell);
  size_per = ref_cell_size_per(ref_cell);

  ref_malloc(sent_c2n, size_per * chunk, REF_GLOB);

  if (ref_mpi_once(ref_mpi)) {
    ref_malloc(elements_to_send, ref_mpi_n(ref_mpi), REF_INT);
    ref_malloc(start_to_send, ref_mpi_n(ref_mpi), REF_INT);
    ref_malloc(c2n, size_per * chunk, REF_GLOB);
    ref_malloc(c2n_int, (node_per + 1) * chunk, REF_INT);
    ref_malloc(c2n_long, (node_per + 1) * chunk, REF_LONG);
    ref_malloc(dest, chunk, REF_INT);

    ncell_read = 0;
    while (ncell_read < ncell) {
      section_size = MIN(chunk, (REF_INT)(ncell - ncell_read));
      nread = (size_t)(section_size * (1 + node_per));
      if (pad) {
        for (cell = 0; cell < section_size; cell++) {
          REF_INT zero;
          REIS(node_per,
               fread(&(c2n_int[(node_per + 1) * cell]), sizeof(REF_INT),
                     (size_t)node_per, file),
               "int c2n pad node");
          REIS(1, fread(&(zero), sizeof(REF_INT), 1, file), "int c2n pad zero");
          REIS(1,
               fread(&(c2n_int[node_per + (node_per + 1) * cell]),
                     sizeof(REF_INT), 1, file),
               "int c2n pad tag");
        }
      } else {
        if (version < 4) {
          REIS(nread, fread(c2n_int, sizeof(REF_INT), nread, file), "int c2n");
          for (cell = 0; cell < section_size; cell++)
            for (node = 0; node < size_per; node++)
              c2n[node + size_per * cell] =
                  (REF_GLOB)c2n_int[node + (node_per + 1) * cell];
        } else {
          REIS(nread, fread(c2n_long, sizeof(REF_LONG), nread, file),
               "long c2n");
          for (cell = 0; cell < section_size; cell++)
            for (node = 0; node < size_per; node++)
              c2n[node + size_per * cell] =
                  (REF_GLOB)c2n_long[node + (node_per + 1) * cell];
        }
      }
      for (cell = 0; cell < section_size; cell++)
        for (node = 0; node < node_per; node++) c2n[node + size_per * cell]--;

      if (REF_CELL_PYR == ref_cell_type(ref_cell)) {
        REF_GLOB n0, n1, n2, n3, n4;
        /* convention: square basis is 0-1-2-3
           (oriented counter clockwise like trias) and top vertex is 4 */
        for (cell = 0; cell < section_size; cell++) {
          n0 = c2n[0 + size_per * cell];
          n1 = c2n[1 + size_per * cell];
          n2 = c2n[2 + size_per * cell];
          n3 = c2n[3 + size_per * cell];
          n4 = c2n[4 + size_per * cell];
          c2n[0 + size_per * cell] = n0;
          c2n[3 + size_per * cell] = n1;
          c2n[4 + size_per * cell] = n2;
          c2n[1 + size_per * cell] = n3;
          c2n[2 + size_per * cell] = n4;
        }
      }

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
    ref_free(c2n_long);
    ref_free(c2n_int);
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

static REF_STATUS ref_part_meshb_cell_bcast(REF_CELL ref_cell, REF_GLOB ncell,
                                            REF_NODE ref_node, REF_INT version,
                                            FILE *file) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_GLOB ncell_read;
  REF_INT chunk;
  REF_INT section_size;
  REF_GLOB *c2n;
  REF_INT *c2n_int;
  REF_LONG *c2n_long;
  REF_INT node_per, size_per;
  REF_INT cell, node, local, new_cell;
  REF_BOOL have_all_nodes, one_node_local;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  size_t nread;

  chunk = (REF_INT)MAX(1000000, ncell / (REF_LONG)ref_mpi_n(ref_mpi));

  size_per = ref_cell_size_per(ref_cell);
  node_per = ref_cell_node_per(ref_cell);

  ref_malloc(c2n, size_per * chunk, REF_GLOB);
  ref_malloc(c2n_int, (1 + node_per) * chunk, REF_INT);
  ref_malloc(c2n_long, (1 + node_per) * chunk, REF_LONG);

  ncell_read = 0;
  while (ncell_read < ncell) {
    section_size = MIN(chunk, (REF_INT)(ncell - ncell_read));
    if (ref_mpi_once(ref_mpi)) {
      nread = (size_t)(section_size * (1 + node_per));
      if (version < 4) {
        REIS(nread, fread(c2n_int, sizeof(REF_INT), nread, file), "int c2n");
        for (cell = 0; cell < section_size; cell++)
          for (node = 0; node < size_per; node++)
            c2n[node + size_per * cell] =
                (REF_GLOB)c2n_int[node + (node_per + 1) * cell];
      } else {
        REIS(nread, fread(c2n_long, sizeof(REF_LONG), nread, file), "long c2n");
        for (cell = 0; cell < section_size; cell++)
          for (node = 0; node < size_per; node++)
            c2n[node + size_per * cell] =
                (REF_GLOB)c2n_long[node + (node_per + 1) * cell];
      }
      for (cell = 0; cell < section_size; cell++)
        for (node = 0; node < node_per; node++) c2n[node + size_per * cell]--;
    }
    RSS(ref_mpi_bcast(ref_mpi, c2n, size_per * section_size, REF_GLOB_TYPE),
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
          nodes[node_per] = (REF_INT)c2n[node_per + size_per * cell];
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

  free(c2n_long);
  free(c2n_int);
  free(c2n);

  return REF_SUCCESS;
}

static REF_STATUS ref_part_meshb(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                                 const char *filename) {
  REF_BOOL verbose = REF_FALSE;
  REF_INT version, dim;
  REF_BOOL available;
  REF_FILEPOS next_position = REF_EMPTY;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_GEOM ref_geom;
  FILE *file;
  REF_BOOL swap_endian = REF_FALSE;
  REF_LONG nnode;
  REF_INT group, keyword_code;
  REF_CELL ref_cell;
  REF_LONG ncell;
  REF_INT type, geom_keyword;
  REF_LONG ngeom;
  REF_INT cad_data_keyword;
  REF_BOOL pad = REF_FALSE;

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
  RSS(ref_mpi_bcast(ref_mpi, &version, 1, REF_INT_TYPE), "bcast");
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
    RSS(ref_part_meshb_long(file, version, &nnode), "nnode");
    if (verbose) printf("nnode %ld\n", nnode);
  }
  RSS(ref_mpi_bcast(ref_mpi, &nnode, 1, REF_LONG_TYPE), "bcast");
  RSS(ref_part_node(file, swap_endian, version, ref_grid_twod(ref_grid),
                    ref_node, nnode),
      "part node");
  if (ref_grid_once(ref_grid))
    REIS(next_position, ftello(file), "vertex file location");

  each_ref_grid_all_ref_cell(ref_grid, group, ref_cell) {
    if (ref_grid_once(ref_grid)) {
      RSS(ref_cell_meshb_keyword(ref_cell, &keyword_code), "kw");
      RSS(ref_import_meshb_jump(file, version, key_pos, keyword_code,
                                &available, &next_position),
          "jump");
      if (available) {
        RSS(ref_part_meshb_long(file, version, &ncell), "ncell");
        if (verbose) printf("group %d ncell %ld\n", group, ncell);
      }
    }
    RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");
    if (available) {
      RSS(ref_mpi_bcast(ref_mpi, &ncell, 1, REF_LONG_TYPE), "bcast");
      RSS(ref_part_meshb_cell(ref_cell, ncell, ref_node, nnode, version, pad,
                              file),
          "part cell");
      if (ref_grid_once(ref_grid))
        REIS(next_position, ftello(file), "cell file location");
    }
  }

  each_ref_type(ref_geom, type) {
    if (ref_grid_once(ref_grid)) {
      geom_keyword = 40 + type;
      RSS(ref_import_meshb_jump(file, version, key_pos, geom_keyword,
                                &available, &next_position),
          "jump");
      if (available) {
        RSS(ref_part_meshb_long(file, version, &ngeom), "ngeom");
        if (verbose) printf("type %d ngeom %ld\n", type, ngeom);
      }
    }
    RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");
    if (available) {
      RSS(ref_mpi_bcast(ref_mpi, &ngeom, 1, REF_LONG_TYPE), "bcast");
      RSS(ref_part_meshb_geom_bcast(ref_geom, ngeom, type, ref_node, version,
                                    file),
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
      RSS(ref_part_meshb_size(file, version,
                              &(ref_geom_cad_data_size(ref_geom))),
          "cad data size");
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
    if (!ref_grid_once(ref_grid)) {
      /* safe non-NULL free, if already allocated, to prevent mem leaks */
      ref_free(ref_geom_cad_data(ref_geom));
      ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                        ref_geom_cad_data_size(ref_geom), REF_BYTE);
    }
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
  REF_INT version = REF_EMPTY, dim = REF_EMPTY;
  REF_BOOL available;
  REF_FILEPOS next_position;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  REF_INT cad_data_keyword;
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

  if (ref_grid_once(ref_grid)) {
    cad_data_keyword = 126; /* GmfByteFlow */
    RSS(ref_import_meshb_jump(file, version, key_pos, cad_data_keyword,
                              &available, &next_position),
        "jump");
    if (available) {
      RSS(ref_part_meshb_size(file, version,
                              &(ref_geom_cad_data_size(ref_geom))),
          "cad data size");
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
    if (!ref_grid_once(ref_grid)) {
      /* safe non-NULL free, if already allocated, to prevent mem leaks */
      ref_free(ref_geom_cad_data(ref_geom));
      ref_malloc_size_t(ref_geom_cad_data(ref_geom),
                        ref_geom_cad_data_size(ref_geom), REF_BYTE);
    }
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
  REF_FILEPOS next_position = REF_EMPTY;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  REF_INT type, geom_keyword;
  REF_LONG ngeom;
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
  RSS(ref_mpi_bcast(ref_mpi, &version, 1, REF_INT_TYPE), "bcast");

  RSS(ref_geom_initialize(ref_geom), "clear out previous assoc");

  each_ref_type(ref_geom, type) {
    if (ref_grid_once(ref_grid)) {
      geom_keyword = 40 + type;
      RSS(ref_import_meshb_jump(file, version, key_pos, geom_keyword,
                                &available, &next_position),
          "jump");
      if (available) {
        RSS(ref_part_meshb_long(file, version, &ngeom), "ngeom");
        if (verbose) printf("type %d ngeom %ld\n", type, ngeom);
      }
    }
    RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");
    if (available) {
      RSS(ref_mpi_bcast(ref_mpi, &ngeom, 1, REF_LONG_TYPE), "bcast");
      RSS(ref_part_meshb_geom_bcast(ref_geom, ngeom, type, ref_node, version,
                                    file),
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
  REF_FILEPOS next_position = REF_EMPTY;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  REF_LONG ncell;
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
  RSS(ref_mpi_bcast(ref_mpi, &version, 1, REF_INT_TYPE), "bcast");

  RSS(ref_cell_free(ref_grid_edg(ref_grid)), "clear out edge");
  RSS(ref_cell_create(&ref_grid_edg(ref_grid), REF_CELL_EDG), "edg");

  if (ref_grid_once(ref_grid)) {
    RSS(ref_import_meshb_jump(file, version, key_pos, 5, &available,
                              &next_position),
        "jump");
    if (available) {
      RSS(ref_part_meshb_long(file, version, &ncell), "ncell");
      if (verbose) printf("nedge %ld\n", ncell);
    }
  }
  RSS(ref_mpi_bcast(ref_mpi, &available, 1, REF_INT_TYPE), "bcast");

  RAS(available, "no edge available in meshb");

  RSS(ref_mpi_bcast(ref_mpi, &ncell, 1, REF_GLOB_TYPE), "bcast");

  RSS(ref_part_meshb_cell_bcast(ref_grid_edg(ref_grid), ncell, ref_node,
                                version, file),
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

  REF_INT version = 0;
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

  RSS(ref_part_node(file, swap_endian, version, REF_FALSE, ref_node, nnode),
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

static REF_STATUS ref_part_avm(REF_GRID *ref_grid_ptr, REF_MPI ref_mpi,
                               const char *filename) {
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_BOOL verbose = REF_TRUE;
  FILE *file;
  REF_INT dim;
  REF_LONG nnode;
  REF_LONG ntet;
  REF_LONG ntri;

  RSS(ref_grid_create(ref_grid_ptr, ref_mpi), "create grid");
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  file = NULL;
  if (ref_mpi_once(ref_mpi)) {
    int i, length, magic, revision, meshes, precision;
    char letter;
    char mesh_type[8];
    char coordinate_system[7];
    double model_scale;
    char units[3];
    double reference[7];
    int refined;
    int nnodes, nfaces, ncells;
    int max_nodes_per_face;
    int max_nodes_per_cell;
    int max_faces_per_cell;
    char element_scheme[33];
    int face_polynomial_order;
    int cell_polynomial_order;
    int boundary_patches;
    int ntet_int, nhex, npri, npyr;
    int ntri_int, ntri2, nqua, nqua2;
    int zeros[5];
    char patch_label[33];
    char patch_type[17];
    int patch_id;

    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");
    length = 6;
    for (i = 0; i < length; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "letter");
      if (verbose) printf("%c", letter);
    }
    if (verbose) printf("\n");
    REIS(1, fread(&magic, sizeof(magic), 1, file), "magic");
    if (verbose) printf("%d magic\n", magic);
    REIS(1, magic, "magic");
    REIS(1, fread(&revision, sizeof(revision), 1, file), "revision");
    if (verbose) printf("%d revision\n", revision);
    REIS(2, revision, "revision");
    REIS(1, fread(&meshes, sizeof(meshes), 1, file), "meshes");
    if (verbose) printf("%d meshes\n", meshes);
    REIS(1, meshes, "meshes");
    length = 128;
    for (i = 0; i < length; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "letter");
      if (verbose) printf("%c", letter);
    }
    if (verbose) printf("\n");
    REIS(1, fread(&precision, sizeof(precision), 1, file), "precision");
    if (verbose) printf("%d precision\n", precision);
    REIS(2, precision, "precision");
    REIS(1, fread(&dim, sizeof(dim), 1, file), "dim");
    if (verbose) printf("%d dim\n", dim);
    RAS(2 <= dim && dim <= 3, "dim");
    REIS(1, fread(&length, sizeof(length), 1, file), "length");
    if (verbose) printf("%d description length\n", length);
    for (i = 0; i < length; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "letter");
      if (verbose) printf("%c", letter);
    }
    if (verbose) printf("\n");
    /* mesh name */
    length = 128;
    for (i = 0; i < length; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "letter");
      if (verbose) printf("%c", letter);
    }
    if (verbose) printf("\n");
    REIS(7, fread(mesh_type, sizeof(char), 7, file), "letter");
    mesh_type[7] = '\0';
    if (verbose) printf("%s", mesh_type);
    REIS(0, strcmp("unstruc", mesh_type), "mesh type");
    length = 128 - 7;
    for (i = 0; i < length; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "letter");
      if (verbose) printf("%c", letter);
    }
    if (verbose) printf("\n");
    /* mesh generator */
    length = 128;
    for (i = 0; i < length; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "letter");
      if (verbose) printf("%c", letter);
    }
    if (verbose) printf("\n");
    REIS(6, fread(coordinate_system, sizeof(char), 6, file),
         "coordinate_system");
    coordinate_system[6] = '\0';
    if (verbose) printf("%s", coordinate_system);
    RSS(ref_grid_parse_coordinate_system(ref_grid, coordinate_system),
        "parse coordinate_system");
    length = 128 - 6;
    for (i = 0; i < length; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "letter");
      if (verbose) printf("%c", letter);
    }
    if (verbose) printf("\n");
    REIS(1, fread(&model_scale, sizeof(model_scale), 1, file), "model_scale");
    if (verbose) printf("%f model_scale\n", model_scale);
    RWDS(1, model_scale, -1, "model_scale");
    REIS(2, fread(units, sizeof(char), 2, file), "units");
    units[2] = '\0';
    if (verbose) printf("%s", units);
    RSS(ref_grid_parse_unit(ref_grid, units), "parse unit");
    length = 128 - 2;
    for (i = 0; i < length; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "letter");
      if (verbose) printf("%c", letter);
    }
    if (verbose) printf("\n");
    REIS(7, fread(reference, sizeof(double), 7, file), "letter");
    for (i = 0; i < 7; i++) {
      ref_grid_reference(ref_grid, i) = reference[i];
      if (verbose) printf("%f reference %d\n", reference[i], i);
    }
    /* reference point description */
    length = 128;
    for (i = 0; i < length; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "letter");
      if (verbose) printf("%c", letter);
    }
    if (verbose) printf("\n");
    REIS(1, fread(&refined, sizeof(refined), 1, file), "refined");
    if (verbose) printf("%d refined\n", refined);
    REIS(0, refined, "refined");
    /* mesh description */
    length = 128;
    for (i = 0; i < length; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "letter");
      if (verbose) printf("%c", letter);
    }
    if (verbose) printf("\n");
    REIS(1, fread(&nnodes, sizeof(nnodes), 1, file), "nnodes");
    REIS(1, fread(&nfaces, sizeof(nfaces), 1, file), "nfaces");
    REIS(1, fread(&ncells, sizeof(ncells), 1, file), "ncells");
    if (verbose)
      printf("%d nnodes %d nfaces %d ncells\n", nnodes, nfaces, ncells);
    REIS(1, fread(&max_nodes_per_face, sizeof(max_nodes_per_face), 1, file),
         "max_nodes_per_face");
    REIS(1, fread(&max_nodes_per_cell, sizeof(max_nodes_per_cell), 1, file),
         "max_nodes_per_cell");
    REIS(1, fread(&max_faces_per_cell, sizeof(max_faces_per_cell), 1, file),
         "max_faces_per_cell");
    if (verbose)
      printf(
          "%d max_nodes_per_face %d max_faces_per_face %d max_faces_per_cell\n",
          max_nodes_per_face, max_nodes_per_cell, max_faces_per_cell);
    REIS(32, fread(element_scheme, sizeof(char), 32, file), "element_scheme");
    element_scheme[32] = '\0';
    if (verbose) printf("%s\n", element_scheme);
    REIS(0, strcmp("uniform", element_scheme), "element_scheme");
    REIS(1,
         fread(&face_polynomial_order, sizeof(face_polynomial_order), 1, file),
         "face_polynomial_order");
    REIS(1,
         fread(&cell_polynomial_order, sizeof(cell_polynomial_order), 1, file),
         "cell_polynomial_order");
    if (verbose)
      printf("%d face_polynomial_order %d cell_polynomial_order\n",
             face_polynomial_order, cell_polynomial_order);
    REIS(1, face_polynomial_order, "face_polynomial_order");
    REIS(1, cell_polynomial_order, "cell_polynomial_order");
    REIS(1, fread(&boundary_patches, sizeof(boundary_patches), 1, file),
         "boundary_patches");
    if (verbose) printf("%d boundary_patches\n", boundary_patches);
    REIS(1, fread(&nhex, sizeof(nhex), 1, file), "nhex");
    REIS(1, fread(&ntet_int, sizeof(ntet_int), 1, file), "ntet_int");
    ntet = (REF_LONG)ntet_int;
    REIS(1, fread(&npri, sizeof(npri), 1, file), "npri");
    REIS(1, fread(&npyr, sizeof(npyr), 1, file), "npyr");
    if (verbose)
      printf("%d nhex %d ntet %d npri %d npyr\n", nhex, ntet_int, npri, npyr);
    REIS(0, nhex, "cant do hex");
    REIS(0, npri, "cant do prism");
    REIS(0, npyr, "cant do pyramid");
    REIS(ncells, ntet_int, "ncells does not match ntet");
    REIS(1, fread(&ntri_int, sizeof(ntri_int), 1, file), "ntri_int");
    ntri = (REF_LONG)ntri_int;
    REIS(1, fread(&ntri2, sizeof(ntri2), 1, file), "ntri2");
    REIS(1, fread(&nqua, sizeof(nqua), 1, file), "nqua");
    REIS(1, fread(&nqua2, sizeof(nqua2), 1, file), "nqua2");
    if (verbose)
      printf("%d ntri %d ntri %d nqua %d nqua\n", ntri_int, ntri2, nqua, nqua2);
    REIS(ntri, ntri2, "ntri mismatch");
    REIS(0, nqua, "cant do quad");
    REIS(nqua, nqua2, "nquad mismatch");
    REIS(5, fread(zeros, sizeof(int), 5, file), "zeros");
    for (i = 0; i < 5; i++) {
      REIS(0, zeros[i], "zeros not zero");
    }
    for (i = 0; i < boundary_patches; i++) {
      REIS(32, fread(patch_label, sizeof(char), 32, file), "patch_label");
      patch_label[32] = '\0';
      REIS(16, fread(patch_type, sizeof(char), 16, file), "patch_type");
      patch_type[16] = '\0';
      REIS(1, fread(&patch_id, sizeof(patch_id), 1, file), "patch_id");
      if (verbose)
        printf("%s %s %d -> %d\n", patch_label, patch_type, patch_id,
               -patch_id);
    }
    nnode = nnodes;
  }
  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &dim, 1, REF_INT_TYPE), "dim");
  ref_grid_twod(ref_grid) = (2 == dim);

  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid),
                    &ref_grid_coordinate_system(ref_grid), 1, REF_INT_TYPE),
      "coordinate_system");
  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &ref_grid_unit(ref_grid), 1,
                    REF_INT_TYPE),
      "unit");
  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &ref_grid_reference(ref_grid, 0), 7,
                    REF_DBL_TYPE),
      "reference");

  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &nnode, 1, REF_LONG_TYPE), "bcast");
  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &ntri, 1, REF_LONG_TYPE), "bcast");
  RSS(ref_mpi_bcast(ref_grid_mpi(ref_grid), &ntet, 1, REF_LONG_TYPE), "bcast");

  {
    REF_BOOL swap_endian = REF_FALSE;
    REF_INT version = 0;
    REF_BOOL twod = REF_FALSE;
    RSS(ref_part_node(file, swap_endian, version, twod, ref_node, nnode),
        "part node");
  }

  if (ref_grid_twod(ref_grid)) {
    REF_CELL ref_cell = ref_grid_edg(ref_grid);
    REF_INT version = 0;
    REF_INT cell;
    REF_BOOL pad = REF_TRUE;
    RSS(ref_part_meshb_cell(ref_cell, ntri, ref_node, nnode, version, pad,
                            file),
        "read edg as tri");
    /* positive face ids */
    each_ref_cell_valid_cell(ref_cell, cell) {
      ref_cell_c2n(ref_cell, ref_cell_id_index(ref_cell), cell) =
          -ref_cell_c2n(ref_cell, ref_cell_id_index(ref_cell), cell);
    }
  } else {
    REF_CELL ref_cell = ref_grid_tri(ref_grid);
    REF_INT version = 0;
    REF_INT cell;
    REF_BOOL pad = REF_FALSE;
    RSS(ref_part_meshb_cell(ref_cell, ntri, ref_node, nnode, version, pad,
                            file),
        "read tri");
    /* positive face ids */
    each_ref_cell_valid_cell(ref_cell, cell) {
      ref_cell_c2n(ref_cell, ref_cell_id_index(ref_cell), cell) =
          -ref_cell_c2n(ref_cell, ref_cell_id_index(ref_cell), cell);
    }
  }

  {
    REF_FILEPOS conn_offset, faceid_offset;
    REF_BOOL swap_endian = REF_FALSE;
    REF_BOOL sixty_four_bit = REF_FALSE;
    conn_offset = 0;
    if (ref_grid_once(ref_grid)) conn_offset = ftello(file);
    faceid_offset = 0;
    RSS(ref_part_bin_ugrid_cell(ref_grid_tet(ref_grid), ntet, ref_node, nnode,
                                file, conn_offset, faceid_offset, swap_endian,
                                sixty_four_bit),
        "read tet");
  }

  if (ref_grid_once(ref_grid)) REIS(0, fclose(file), "close file");
  return REF_SUCCESS;
}

static REF_STATUS ref_part_metric_solb(REF_NODE ref_node,
                                       const char *filename) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_FILEPOS next_position = REF_EMPTY;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  FILE *file;
  REF_INT chunk;
  REF_DBL *metric;
  REF_LONG nnode_read;
  REF_INT section_size;
  REF_INT node, local;
  REF_BOOL available;
  REF_INT version, dim, ntype, type;
  REF_GLOB global;
  REF_LONG nnode;

  file = NULL;
  if (ref_mpi_once(ref_mpi)) {
    RSS(ref_import_meshb_header(filename, &version, key_pos), "header");
    RAS(2 <= version && version <= 4, "unsupported version");
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");
    RSS(ref_import_meshb_jump(file, version, key_pos, 3, &available,
                              &next_position),
        "jump");
    RAS(available, "solb missing dimension");
    REIS(1, fread((unsigned char *)&dim, 4, 1, file), "dim");
    RAS(2 <= dim && dim <= 3, "unsupported dimension");
    RSS(ref_import_meshb_jump(file, version, key_pos, 62, &available,
                              &next_position),
        "jump");
    RAS(available, "SolAtVertices missing");
    RSS(ref_part_meshb_long(file, version, &nnode), "nnode");
    REIS(1, fread((unsigned char *)&ntype, 4, 1, file), "ntype");
    REIS(1, fread((unsigned char *)&type, 4, 1, file), "type");
    REIS(1, ntype, "number of solutions");
    REIS(3, type, "metric solution type");
  }
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), &version, 1, REF_INT_TYPE),
      "bcast version");
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), &dim, 1, REF_INT_TYPE),
      "bcast dim");
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), &nnode, 1, REF_GLOB_TYPE),
      "bcast nnode");

  if ((nnode != ref_node_n_global(ref_node)) &&
      (nnode / 2 != ref_node_n_global(ref_node))) {
    if (ref_mpi_once(ref_mpi))
      printf("file %ld ref_node " REF_GLOB_FMT " %s\n", nnode,
             ref_node_n_global(ref_node), filename);
    THROW("global count mismatch");
  }

  chunk =
      (REF_INT)MAX(100000, nnode / (REF_LONG)ref_mpi_n(ref_node_mpi(ref_node)));
  chunk = (REF_INT)MIN((REF_LONG)chunk, nnode);

  ref_malloc_init(metric, 6 * chunk, REF_DBL, -1.0);

  nnode_read = 0;
  while (nnode_read < nnode) {
    section_size = MIN(chunk, (REF_INT)(nnode - nnode_read));
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
          REIS(1, fread(&(metric[1 + 6 * node]), sizeof(REF_DBL), 1, file),
               "m12");
          REIS(1, fread(&(metric[3 + 6 * node]), sizeof(REF_DBL), 1, file),
               "m22");
          metric[2 + 6 * node] = 0.0; /* m13 */
          metric[4 + 6 * node] = 0.0; /* m23 */
          metric[5 + 6 * node] = 1.0; /* m33 */
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
      global = (REF_LONG)node + nnode_read;
      RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
      if (REF_EMPTY != local) {
        RSS(ref_node_metric_set(ref_node, local, &(metric[6 * node])),
            "set local node met");
      }
      if (2 == dim) {
        global = nnode + (REF_GLOB)node + nnode_read;
        RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
        if (REF_EMPTY != local) {
          RSS(ref_node_metric_set(ref_node, local, &(metric[6 * node])),
              "set local node met");
        }
      }
    }
    nnode_read += (REF_LONG)section_size;
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
  REF_INT dim = REF_EMPTY;

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
      dim = REF_EMPTY;
      while (!feof(file)) {
        status = fscanf(file, "%s", line);
        if (EOF == status) break;
        REIS(1, status, "line read failed");

        if (0 == strcmp("Dimension", line)) {
          REIS(1, fscanf(file, "%d", &dim), "read dim");
        }

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
      RUS(REF_EMPTY, dim, "Dimension keyword missing from .sol metric");
      RAS(found_keyword, "SolAtVertices keyword missing from .sol metric");
    } else {
      nnode = (REF_INT)ref_node_n_global(ref_node);
    }
  }
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), &nnode, 1, REF_INT_TYPE),
      "bcast nnode");

  chunk = (REF_INT)MAX(100000, nnode / ref_mpi_n(ref_node_mpi(ref_node)));
  chunk = (REF_INT)MIN((REF_GLOB)chunk, nnode);

  ref_malloc_init(metric, 6 * chunk, REF_DBL, -1.0);

  nnode_read = 0;
  while (nnode_read < nnode) {
    section_size =
        (REF_INT)MIN((REF_GLOB)chunk, ref_node_n_global(ref_node) - nnode_read);
    if (ref_mpi_once(ref_node_mpi(ref_node))) {
      for (node = 0; node < section_size; node++)
        if (sol_format) {
          if (3 == dim) {
            REIS(6,
                 fscanf(file, "%lf %lf %lf %lf %lf %lf",
                        &(metric[0 + 6 * node]), &(metric[1 + 6 * node]),
                        &(metric[3 + 6 * node]), /* transposed 3,2 */
                        &(metric[2 + 6 * node]), &(metric[4 + 6 * node]),
                        &(metric[5 + 6 * node])),
                 "metric read error");
          } else {
            REIS(3,
                 fscanf(file, "%lf %lf %lf", &(metric[0 + 6 * node]),
                        &(metric[1 + 6 * node]), &(metric[3 + 6 * node])),
                 "metric read error");
            metric[2 + 6 * node] = 0.0; /* m13 */
            metric[4 + 6 * node] = 0.0; /* m23 */
            metric[5 + 6 * node] = 1.0; /* m33 */
          }
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

static REF_STATUS ref_part_plt_string(FILE *file, char *string, int maxlen) {
  int i, letter;
  for (i = 0; i < maxlen; i++) {
    REIS(1, fread(&letter, sizeof(int), 1, file), "plt string letter");
    string[i] = (char)letter;
    if (0 == letter) {
      return REF_SUCCESS;
    }
  }
  return REF_FAILURE;
}

static REF_STATUS ref_part_plt_header(FILE *file, REF_INT *nvar,
                                      REF_LIST zone_type, REF_LIST zone_packing,
                                      REF_LIST zone_nnode,
                                      REF_LIST zone_nelem) {
  char header[9];
  int endian, filetype;
  char title[1024], varname[1024], zonename[1024];
  int var, numvar;
  float zonemarker;
  int parent, strand, notused, zonetype, packing, location, neighbor;
  double solutiontime;
  int miscellaneous, i;
  int numpts, numelem;
  int dim, aux;
  REF_BOOL verbose = REF_FALSE;

  RAS(header == fgets(header, 6, file), "header error");
  header[5] = '\0';
  REIS(0, strncmp(header, "#!TDV", 5), "header '#!TDV' missing")
  RAS(header == fgets(header, 4, file), "version error");
  header[4] = '\0';
  REIS(0, strncmp(header, "112", 5), "expected version '112'")

  REIS(1, fread(&endian, sizeof(int), 1, file), "magic");
  REIS(1, endian, "expected little endian plt");
  REIS(1, fread(&filetype, sizeof(int), 1, file), "filetype");
  REIS(0, filetype, "expected full filetype");

  RSS(ref_part_plt_string(file, title, 1024), "read title");
  if (verbose) printf("plt title '%s'\n", title);

  REIS(1, fread(&numvar, sizeof(int), 1, file), "numvar");
  if (verbose) printf("plt number of variables %d\n", numvar);
  *nvar = numvar;

  for (var = 0; var < numvar; var++) {
    RSS(ref_part_plt_string(file, varname, 1024), "read variable name");
    if (verbose) printf("plt variable name %d '%s'\n", var, varname);
  }

  REIS(1, fread(&zonemarker, sizeof(float), 1, file), "zonemarker");
  while (ABS(299.0 - (double)zonemarker) < 1.0e-7) {
    RSS(ref_part_plt_string(file, zonename, 1024), "read zonename");
    if (verbose) printf("plt zonename '%s'\n", zonename);

    REIS(1, fread(&parent, sizeof(int), 1, file), "parent");
    REIS(1, fread(&strand, sizeof(int), 1, file), "strand");
    REIS(1, fread(&solutiontime, sizeof(double), 1, file), "solutiontime");
    REIS(1, fread(&notused, sizeof(int), 1, file), "notused");
    REIS(-1, notused, "not unused should be -1 plt");
    REIS(1, fread(&zonetype, sizeof(int), 1, file), "zonetype");
    RSS(ref_list_push(zone_type, zonetype), "save zonetype");
    REIS(1, fread(&packing, sizeof(int), 1, file), "packing");
    RSS(ref_list_push(zone_packing, packing), "save packing");
    REIS(1, fread(&location, sizeof(int), 1, file), "location");
    REIS(0, location, "only node data location plt implemented");
    REIS(1, fread(&neighbor, sizeof(int), 1, file), "neighbor");
    REIS(0, neighbor, "no face neighbor  plt implemented");

    /* there seem to be numvar extra zeros in USM3D volume files */
    /* assume numpts is first nonzero */
    miscellaneous = 0;
    for (i = 0; (i < (numvar + 1)) && (0 == miscellaneous); i++) {
      REIS(1, fread(&miscellaneous, sizeof(int), 1, file), "mystery");
    }
    numpts = miscellaneous;

    RSS(ref_list_push(zone_nnode, numpts), "save nnode");
    REIS(1, fread(&numelem, sizeof(int), 1, file), "numelem");
    RSS(ref_list_push(zone_nelem, numelem), "save nelem");

    for (i = 0; i < 3; i++) {
      REIS(1, fread(&dim, sizeof(int), 1, file), "dim");
      REIS(0, dim, "dim nonzero plt");
    }
    REIS(1, fread(&aux, sizeof(int), 1, file), "aux");
    REIS(0, aux, "aux nonzero plt");

    REIS(1, fread(&zonemarker, sizeof(float), 1, file), "zonemarker");
  }

  RWDS(357.0, (double)zonemarker, -1.0, "end of header marker expected");

  return REF_SUCCESS;
}

static REF_STATUS ref_part_plt_data(FILE *file, REF_INT nvar,
                                    REF_LIST zone_type, REF_LIST zone_packing,
                                    REF_LIST zone_nnode, REF_LIST zone_nelem,
                                    REF_INT *length, REF_DBL **soln) {
  float zonemarker;
  int dataformat;
  REF_INT i, node, elem, node_per;
  int passive, sharing, conn, c2n;
  double minval, maxval;
  REF_INT zonetype, packing, nnode, nelem;
  float var_float;
  double var_double;
  REF_LIST dataformats = NULL;
  REF_BOOL verbose = REF_FALSE;

  RSS(ref_list_create(&dataformats), "dataformats");

  RSS(ref_list_shift(zone_type, &zonetype), "zonetype");
  RSS(ref_list_shift(zone_packing, &packing), "zone packing");
  RSS(ref_list_shift(zone_nnode, &nnode), "zone node size");
  RSS(ref_list_shift(zone_nelem, &nelem), "zone elem size");

  switch (zonetype) {
    case 1: /* FELINESEG */
      node_per = 2;
      break;
    case 2: /* FETRIANGLE */
      node_per = 3;
      break;
    case 3: /* FEQUADRILATERAL */
      node_per = 4;
      break;
    case 4: /* FETETRAHEDRON */
      node_per = 4;
      break;
    case 5: /* FEBRICK */
      node_per = 8;
      break;
    default:
      printf("zonetype = %d\n", zonetype);
      THROW("unknown tecplot plt zonetype read");
  }
  if (verbose)
    printf("zone type %d node_per %d packing %d nnode %d nelem %d\n", zonetype,
           node_per, packing, nnode, nelem);

  *length = nnode;

  REIS(1, fread(&zonemarker, sizeof(float), 1, file), "zonemarker");
  RWDS(299.0, (double)zonemarker, -1.0, "start of data header expected");

  for (i = 0; i < nvar; i++) {
    REIS(1, fread(&dataformat, sizeof(int), 1, file), "dim");
    RSS(ref_list_push(dataformats, dataformat), "save dataformat");
  }
  REIS(1, fread(&passive, sizeof(int), 1, file), "dim");
  if (1 == passive) {
    for (i = 0; i < nvar; i++) {
      REIS(1, fread(&passive, sizeof(int), 1, file), "dim");
      REIS(0, passive, "passive variable nonzero plt");
    }
  }
  REIS(1, fread(&sharing, sizeof(int), 1, file), "dim");
  if (1 == sharing) {
    for (i = 0; i < nvar; i++) {
      REIS(1, fread(&sharing, sizeof(int), 1, file), "dim");
      REIS(-1, sharing, "variable sharing not implemented plt");
    }
  }
  REIS(1, fread(&conn, sizeof(int), 1, file), "dim");
  REIS(-1, conn, "connectivity sharing not implemented plt");

  for (i = 0; i < nvar; i++) {
    REIS(1, fread(&minval, sizeof(double), 1, file), "dim");
    REIS(1, fread(&maxval, sizeof(double), 1, file), "dim");
  }

  ref_malloc(*soln, nvar * nnode, REF_DBL);

  /* Data packing. 0 = Block 1 = Point */
  if (0 == packing) {
    for (i = 0; i < nvar; i++) {
      for (node = 0; node < nnode; node++) {
        dataformat = ref_list_value(dataformats, i);
        switch (dataformat) {
          case 1: /* float */
            REIS(1, fread(&var_float, sizeof(float), 1, file), "float");
            (*soln)[i + nvar * node] = (double)var_float;
            break;
          case 2: /* double */
            REIS(1, fread(&var_double, sizeof(double), 1, file), "double");
            (*soln)[i + nvar * node] = var_double;
            break;
          default:
            printf("dataformat = %d\n", dataformat);
            RSS(REF_IMPLEMENT, "implement tecplot plt dataformat read");
        }
      }
    }
  } else {
    for (node = 0; node < nnode; node++) {
      for (i = 0; i < nvar; i++) {
        dataformat = ref_list_value(dataformats, i);
        switch (dataformat) {
          case 1: /* float */
            REIS(1, fread(&var_float, sizeof(float), 1, file), "float");
            (*soln)[i + nvar * node] = (double)var_float;
            break;
          case 2: /* double */
            REIS(1, fread(&var_double, sizeof(double), 1, file), "double");
            (*soln)[i + nvar * node] = var_double;
            break;
          default:
            printf("dataformat = %d\n", dataformat);
            RSS(REF_IMPLEMENT, "implement tecplot plt dataformat read");
        }
      }
    }
  }

  for (elem = 0; elem < nelem; elem++) {
    for (i = 0; i < node_per; i++) {
      REIS(1, fread(&c2n, sizeof(int), 1, file), "dim");
    }
  }

  ref_list_free(dataformats);

  return REF_SUCCESS;
}

static REF_STATUS ref_part_scalar_plt(REF_GRID ref_grid, REF_INT *ldim,
                                      REF_DBL **scalar, const char *filename) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  FILE *file = NULL;
  REF_INT nvar;
  REF_LIST zone_type = NULL, zone_packing = NULL, zone_nnode = NULL,
           zone_nelem = NULL, touching;
  REF_INT zone, nzone, length, node, i, point;
  REF_DBL *soln = NULL;
  REF_SEARCH ref_search;
  REF_DBL radius, position[3], dist, best_dist;
  REF_INT best, item;
  REF_BOOL verbose = REF_FALSE;

  RSS(ref_search_create(&ref_search, ref_node_n(ref_node)), "create search");
  each_ref_node_valid_node(ref_node, node) {
    radius = 0.0;
    RSS(ref_search_insert(ref_search, node, ref_node_xyz_ptr(ref_node, node),
                          radius),
        "ins");
  }

  if (ref_mpi_once(ref_mpi)) {
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    RSS(ref_list_create(&zone_type), "zonetype list");
    RSS(ref_list_create(&zone_packing), "packing list");
    RSS(ref_list_create(&zone_nnode), "nnode list");
    RSS(ref_list_create(&zone_nelem), "nelem list");

    RSS(ref_part_plt_header(file, &nvar, zone_type, zone_packing, zone_nnode,
                            zone_nelem),
        "parse header");
    nzone = ref_list_n(zone_nnode);

    {
      REF_BOOL force_block = REF_FALSE;
      each_ref_list_item(zone_packing, item) {
        force_block = force_block || (1 == ref_list_value(zone_packing, item));
        ref_list_value(zone_packing, item) = 0;
      }
      if (force_block) printf("data packing set to block, was point\n");
    }
  }
  RSS(ref_mpi_bcast(ref_mpi, &nvar, 1, REF_INT_TYPE), "b nvar");
  RSS(ref_mpi_bcast(ref_mpi, &nzone, 1, REF_INT_TYPE), "b nzone");

  *ldim = nvar - 3;
  ref_malloc_init(*scalar, (*ldim) * ref_node_max(ref_node), REF_DBL, -999.0);

  RSS(ref_list_create(&touching), "touching list");
  for (zone = 0; zone < nzone; zone++) {
    if (ref_mpi_once(ref_mpi)) {
      RSS(ref_part_plt_data(file, nvar, zone_type, zone_packing, zone_nnode,
                            zone_nelem, &length, &soln),
          "read data");
      RSS(ref_mpi_bcast(ref_mpi, &length, 1, REF_INT_TYPE), "b length");
      if (length > 0) {
        RSS(ref_mpi_bcast(ref_mpi, soln, nvar * length, REF_DBL_TYPE),
            "b soln");
      }
    } else {
      RSS(ref_mpi_bcast(ref_mpi, &length, 1, REF_INT_TYPE), "b length");
      if (length > 0) {
        ref_malloc(soln, nvar * length, REF_DBL);
        RSS(ref_mpi_bcast(ref_mpi, soln, nvar * length, REF_DBL_TYPE),
            "b soln");
      } else {
        soln = NULL;
      }
    }
    for (point = 0; point < length; point++) {
      if (verbose) {
        printf("point %d", point);
        for (i = 0; i < nvar; i++) {
          printf(" %e", soln[i + nvar * point]);
        }
        printf("\n");
      }

      if (ref_grid_twod(ref_grid)) {
        position[0] = soln[0 + nvar * point];
        position[1] = soln[2 + nvar * point];
        position[2] = 0.0;
      } else {
        for (i = 0; i < 3; i++) {
          position[i] = soln[i + nvar * point];
        }
      }
      /* single precision */
      radius =
          100.0 * 1.0e-8 *
          sqrt(pow(position[0], 2) + pow(position[1], 2) + pow(position[2], 2));
      RSS(ref_search_touching(ref_search, touching, position, radius),
          "search tree");
      best_dist = 1.0e+200;
      best = REF_EMPTY;
      each_ref_list_item(touching, item) {
        node = ref_list_value(touching, item);
        dist = sqrt(pow(ref_node_xyz(ref_node, 0, node) - position[0], 2) +
                    pow(ref_node_xyz(ref_node, 1, node) - position[1], 2) +
                    pow(ref_node_xyz(ref_node, 2, node) - position[2], 2));
        if (dist < best_dist) {
          best_dist = dist;
          best = node;
        }
      }
      if (verbose) printf("best %d %e\n", best, best_dist);
      if (REF_EMPTY != best) {
        for (i = 3; i < nvar; i++) {
          (*scalar)[(i - 3) + (*ldim) * best] = soln[i + nvar * point];
        }
      }
      RSS(ref_list_erase(touching), "erase");
    }
    ref_free(soln);
  }
  RSS(ref_list_free(touching), "free touching");

  if (ref_mpi_once(ref_mpi)) {
    fclose(file);
    ref_list_free(zone_nelem);
    ref_list_free(zone_nnode);
    ref_list_free(zone_packing);
    ref_list_free(zone_type);
  }

  RSS(ref_search_free(ref_search), "free search");

  return REF_SUCCESS;
}

static REF_STATUS ref_part_scalar_cell_restart_sol(REF_GRID ref_grid,
                                                   REF_INT *ldim,
                                                   REF_DBL **scalar,
                                                   const char *filename) {
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  FILE *file;
  REF_INT ncell, i, cell_node;
  REF_INT node, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT *hits;

  RAS(!ref_mpi_para(ref_mpi), "only implemented for single core");

  file = fopen(filename, "r");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  *ldim = 5;
  REIS(1, fscanf(file, "%d", &ncell), "read ncell");
  if (ncell != ref_cell_n(ref_cell)) {
    printf("file %d ref_cell %d\n", ncell, ref_cell_n(ref_cell));
    THROW("ERROR: global count mismatch");
  }

  ref_malloc_init(*scalar, (*ldim) * ref_node_max(ref_node), REF_DBL, 0.0);
  ref_malloc_init(hits, ref_node_max(ref_node), REF_INT, 0.0);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    REF_DBL rho, u, v, p;
    REIS(4, fscanf(file, "%lf %lf %lf %lf", &rho, &u, &v, &p),
         "cell data read");
    each_ref_cell_cell_node(ref_cell, cell_node) {
      (*scalar)[0 + (*ldim) * nodes[cell_node]] += rho;
      (*scalar)[1 + (*ldim) * nodes[cell_node]] += u;
      (*scalar)[3 + (*ldim) * nodes[cell_node]] += v;
      (*scalar)[4 + (*ldim) * nodes[cell_node]] += p;
      hits[nodes[cell_node]] += 1;
    }
  }

  each_ref_node_valid_node(ref_node, node) {
    RAS(hits[node] > 0, "zero hits for node");
    for (i = 0; i < 5; i++) {
      (*scalar)[i + (*ldim) * node] /= (REF_DBL)hits[node];
    }
  }
  ref_free(hits);

  REIS(0, fclose(file), "close file");

  return REF_SUCCESS;
}

static REF_STATUS ref_part_scalar_rst(REF_NODE ref_node, REF_INT *ldim,
                                      REF_DBL **scalar, const char *filename) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  FILE *file;
  REF_BOOL verbose = REF_TRUE;
  int dim, variables, step, steps, dof;
  REF_INT chunk;
  REF_DBL *data;
  REF_INT section_size;
  REF_GLOB global;
  REF_INT node, local;
  REF_LONG nnode_read;

  file = NULL;
  if (ref_mpi_once(ref_mpi)) {
    int i, length, version, doubles;
    char letter;
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");
    REIS(1, fread(&length, sizeof(length), 1, file), "dim");
    if (verbose) printf("%d length of magic string\n", length);
    REIS(8, length, "magic string length");
    for (i = 0; i < length; i++) {
      REIS(1, fread(&letter, sizeof(letter), 1, file), "dim");
      if (verbose) printf("%c", letter);
    }
    if (verbose) printf("\n");
    REIS(1, fread(&version, sizeof(version), 1, file), "dim");
    if (verbose) printf("%d version\n", version);
    REIS(2, version, "version");
    REIS(1, fread(&dim, sizeof(dim), 1, file), "dim");
    if (verbose) printf("%d dim\n", dim);
    RAS(2 <= dim && dim <= 3, "dim");
    REIS(1, fread(&variables, sizeof(variables), 1, file), "variables");
    if (verbose) printf("%d variables\n", variables);
    REIS(1, fread(&steps, sizeof(steps), 1, file), "steps");
    if (verbose) printf("%d steps\n", steps);
    REIS(1, fread(&dof, sizeof(dof), 1, file), "dof");
    if (verbose) printf("%d dof\n", dof);
    REIS(1, fread(&doubles, sizeof(doubles), 1, file), "doubles");
    if (verbose) printf("%d doubles\n", doubles);
    REIS(0, doubles, "expected zero doubles");
    /* assert zero doubles, skip misc metadata (timestep) */
  }
  RSS(ref_mpi_bcast(ref_mpi, &dim, 1, REF_INT_TYPE), "bcast dim");
  RSS(ref_mpi_bcast(ref_mpi, &variables, 1, REF_INT_TYPE), "bcast dim");
  RSS(ref_mpi_bcast(ref_mpi, &steps, 1, REF_INT_TYPE), "bcast dim");
  RSS(ref_mpi_bcast(ref_mpi, &dof, 1, REF_INT_TYPE), "bcast dim");

  if (ref_node_n_global(ref_node) != dof) {
    if (ref_mpi_once(ref_mpi)) {
      printf("file %d ref_node " REF_GLOB_FMT " %s\n", dof,
             ref_node_n_global(ref_node), filename);
    }
    THROW("ERROR: global count mismatch");
  }

  *ldim = variables * steps;
  ref_malloc(*scalar, (*ldim) * ref_node_max(ref_node), REF_DBL);

  chunk =
      (REF_INT)MAX(100000, dof / (REF_LONG)ref_mpi_n(ref_node_mpi(ref_node)));
  chunk = (REF_INT)MIN((REF_LONG)chunk, dof);

  ref_malloc_init(data, variables * chunk, REF_DBL, -1.0);

  for (step = 0; step < steps; step++) {
    int i;
    nnode_read = 0;
    while (nnode_read < dof) {
      section_size = MIN(chunk, (REF_INT)(dof - nnode_read));
      if (ref_mpi_once(ref_node_mpi(ref_node))) {
        REIS((variables)*section_size,
             fread(data, sizeof(REF_DBL), (size_t)(variables * section_size),
                   file),
             "dat");
        RSS(ref_mpi_bcast(ref_node_mpi(ref_node), data, variables * chunk,
                          REF_DBL_TYPE),
            "bcast");
      } else {
        RSS(ref_mpi_bcast(ref_node_mpi(ref_node), data, variables * chunk,
                          REF_DBL_TYPE),
            "bcast");
      }
      for (node = 0; node < section_size; node++) {
        global = node + nnode_read;
        RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
        if (REF_EMPTY != local) {
          for (i = 0; i < variables; i++) {
            (*scalar)[i + step * variables + local * (*ldim)] =
                data[i + node * variables];
          }
        }
      }
      nnode_read += (REF_LONG)section_size;
    }
  }
  ref_free(data);

  if (ref_mpi_once(ref_mpi)) {
    fclose(file);
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_part_scalar_sol(REF_NODE ref_node, REF_INT *ldim,
                                      REF_DBL **scalar, const char *filename) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  FILE *file;
  REF_INT chunk;
  REF_DBL *data;
  REF_INT section_size;
  REF_GLOB global;
  REF_INT node, local;
  REF_INT dim, ntype, type, i;
  REF_LONG nnode = REF_EMPTY, nnode_read;
  REF_INT status;
  char line[1024];
  REF_BOOL found_keyword = REF_FALSE;

  file = NULL;
  if (ref_mpi_once(ref_node_mpi(ref_node))) {
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    dim = REF_EMPTY;
    while (!feof(file)) {
      status = fscanf(file, "%s", line);
      if (EOF == status) break;
      REIS(1, status, "line read failed");

      if (0 == strcmp("Dimension", line)) {
        REIS(1, fscanf(file, "%d", &dim), "read dim");
        printf("dim %d\n", dim);
      }

      if (0 == strcmp("SolAtVertices", line)) {
        REIS(1, fscanf(file, "%ld", &nnode), "read nnode");
        printf("nnode %d\n", (REF_INT)nnode);
        REIS(1, fscanf(file, "%d", &ntype), "read ldim");
        printf("ntype %d\n", ntype);
        *ldim = 0;
        for (i = 0; i < ntype; i++) {
          REIS(1, fscanf(file, "%d", &type), "read type");
          printf("%d type %d\n", i, type);
          if (1 == type) {
            (*ldim) += 1;
          } else if (2 == type) {
            (*ldim) += dim;
          } else {
            printf("item %d type %d\n", i, type);
            THROW("solb type");
          }
        }
        RAS(0 <= fscanf(file, "%*[^1234567890-+.]"), "skip blank line");
        found_keyword = REF_TRUE;
        break;
      }
    }
    RUS(REF_EMPTY, dim, "Dimension keyword missing from .sol metric");
    RAS(found_keyword, "SolAtVertices keyword missing from .sol metric");
    printf("nnode %d ldim %d dim %d\n", (REF_INT)nnode, *ldim, dim);
  }
  RSS(ref_mpi_bcast(ref_mpi, &nnode, 1, REF_GLOB_TYPE), "bcast nnode");
  RSS(ref_mpi_bcast(ref_mpi, ldim, 1, REF_INT_TYPE), "bcast ldim");
  RSS(ref_mpi_bcast(ref_mpi, &dim, 1, REF_INT_TYPE), "bcast dim");

  if ((nnode != ref_node_n_global(ref_node)) &&
      (nnode / 2 != ref_node_n_global(ref_node))) {
    if (ref_mpi_once(ref_mpi))
      printf("file %ld ref_node " REF_GLOB_FMT " %s\n", nnode,
             ref_node_n_global(ref_node), filename);
    if (nnode > ref_node_n_global(ref_node)) {
      if (ref_mpi_once(ref_mpi))
        REF_WHERE("WARNING: global count mismatch, too many");
    } else {
      THROW("ERROR: global count mismatch, too few");
    }
  }

  ref_malloc(*scalar, (*ldim) * ref_node_max(ref_node), REF_DBL);

  chunk =
      (REF_INT)MAX(100000, nnode / (REF_LONG)ref_mpi_n(ref_node_mpi(ref_node)));
  chunk = (REF_INT)MIN((REF_LONG)chunk, nnode);

  ref_malloc_init(data, (*ldim) * chunk, REF_DBL, -1.0);

  nnode_read = 0;
  while (nnode_read < nnode) {
    section_size = MIN(chunk, (REF_INT)(nnode - nnode_read));
    if (ref_mpi_once(ref_node_mpi(ref_node))) {
      for (i = 0; i < ((*ldim) * section_size); i++)
        REIS(1, fscanf(file, "%lf", &(data[i])), "vertex data read error");
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
      if (2 == dim) {
        global = nnode + (REF_GLOB)node + nnode_read;
        RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
        if (REF_EMPTY != local) {
          for (i = 0; i < *ldim; i++) {
            (*scalar)[i + local * (*ldim)] = data[i + node * (*ldim)];
          }
        }
      }
    }
    nnode_read += (REF_LONG)section_size;
  }

  ref_free(data);

  if (ref_mpi_once(ref_node_mpi(ref_node))) {
    REIS(0, fclose(file), "close file");
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_part_scalar_solb(REF_NODE ref_node, REF_INT *ldim,
                                       REF_DBL **scalar, const char *filename) {
  REF_MPI ref_mpi = ref_node_mpi(ref_node);
  REF_FILEPOS next_position = REF_EMPTY;
  REF_FILEPOS key_pos[REF_IMPORT_MESHB_LAST_KEYWORD];
  FILE *file;
  REF_INT chunk;
  REF_DBL *data;
  REF_INT section_size;
  REF_GLOB global;
  REF_INT node, local;
  REF_BOOL available;
  REF_INT version, dim, ntype, type, i;
  REF_LONG nnode, nnode_read;

  file = NULL;
  if (ref_mpi_once(ref_node_mpi(ref_node))) {
    RSS(ref_import_meshb_header(filename, &version, key_pos), "head");
    RAS(2 <= version && version <= 4, "unsupported version");
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");
    RSS(ref_import_meshb_jump(file, version, key_pos, 3, &available,
                              &next_position),
        "jump");
    RAS(available, "solb missing dimension");
    REIS(1, fread((unsigned char *)&dim, 4, 1, file), "dim");
    RAS(2 <= dim && dim <= 3, "unsupported dimension");

    RSS(ref_import_meshb_jump(file, version, key_pos, 62, &available,
                              &next_position),
        "jmp");
    RAS(available, "SolAtVertices missing");
    RSS(ref_part_meshb_long(file, version, &nnode), "nnode");
    REIS(1, fread((unsigned char *)&ntype, 4, 1, file), "ntype");
    *ldim = 0;
    for (i = 0; i < ntype; i++) {
      REIS(1, fread((unsigned char *)&type, 4, 1, file), "type");
      RAB(1 <= type && type <= 2,
          "only types 1 (scalar) or 2 (vector) supported",
          { printf(" %d type\n", type); });
      if (1 == type) (*ldim) += 1;
      if (2 == type) (*ldim) += dim;
    }
  }
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), &version, 1, REF_INT_TYPE),
      "bcast version");
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), &dim, 1, REF_INT_TYPE),
      "bcast dim");
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), &nnode, 1, REF_GLOB_TYPE),
      "bcast nnode");
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), ldim, 1, REF_INT_TYPE),
      "bcast ldim");

  if ((nnode != ref_node_n_global(ref_node)) &&
      (nnode / 2 != ref_node_n_global(ref_node))) {
    if (ref_mpi_once(ref_mpi))
      printf("file %ld ref_node " REF_GLOB_FMT " %s\n", nnode,
             ref_node_n_global(ref_node), filename);
    if (nnode > ref_node_n_global(ref_node)) {
      if (ref_mpi_once(ref_mpi))
        REF_WHERE("WARNING: global count mismatch, too many");
    } else {
      THROW("ERROR: global count mismatch, too few");
    }
  }

  ref_malloc(*scalar, (*ldim) * ref_node_max(ref_node), REF_DBL);

  chunk =
      (REF_INT)MAX(100000, nnode / (REF_LONG)ref_mpi_n(ref_node_mpi(ref_node)));
  chunk = (REF_INT)MIN((REF_LONG)chunk, nnode);

  ref_malloc_init(data, (*ldim) * chunk, REF_DBL, -1.0);

  nnode_read = 0;
  while (nnode_read < nnode) {
    section_size = MIN(chunk, (REF_INT)(nnode - nnode_read));
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
      if (2 == dim) {
        global = nnode + (REF_GLOB)node + nnode_read;
        RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
        if (REF_EMPTY != local) {
          for (i = 0; i < *ldim; i++) {
            (*scalar)[i + local * (*ldim)] = data[i + node * (*ldim)];
          }
        }
      }
    }
    nnode_read += (REF_LONG)section_size;
  }

  ref_free(data);

  if (ref_mpi_once(ref_node_mpi(ref_node))) {
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
  REF_INT version;
  uint64_t version_number;
  uint64_t number_of_fields;
  uint64_t number_of_chars;
  char letter;
  unsigned long field_length, uint_nnode;
  int association;
  REF_INT field;
  REF_BOOL verbose = REF_FALSE;

  nnode = -1;
  next_position = -1;
  file = NULL;
  if (ref_mpi_once(ref_node_mpi(ref_node))) {
    file = fopen(filename, "r");
    if (NULL == (void *)file) printf("unable to open %s\n", filename);
    RNS(file, "unable to open file");

    end_of_string = strlen(filename);
    REIS(0, strcmp(&filename[end_of_string - 5], ".snap"),
         "snap extension expected");

    if (verbose) printf("read snap %s\n", filename);

    REIS(1, fread(&version_number, sizeof(version_number), 1, file),
         "version number");
    RAS(2 <= version_number && version_number <= 3,
        "only versions 2 and 3 supported");
    version = (REF_INT)version_number;

    REIS(1, fread(&number_of_fields, sizeof(number_of_fields), 1, file),
         "number");
    *ldim = (REF_INT)number_of_fields;
  }
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), &version, 1, REF_INT_TYPE),
      "bcast ldim");
  RSS(ref_mpi_bcast(ref_node_mpi(ref_node), ldim, 1, REF_INT_TYPE),
      "bcast ldim");
  ref_malloc(*scalar, (*ldim) * ref_node_max(ref_node), REF_DBL);

  for (field = 0; field < (*ldim); field++) {
    if (ref_mpi_once(ref_node_mpi(ref_node))) {
      if (2 == version) {
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
      } else {
        uint64_t length_remaining;
        uint64_t n_global;
        uint64_t entry_length_uint64;
        uint64_t pair;
        uint64_t count;
        REIS(1, fread(&length_remaining, sizeof(length_remaining), 1, file),
             "length_remaining");
        next_position = (REF_FILEPOS)length_remaining + ftello(file);
        REIS(1, fread(&n_global, sizeof(n_global), 1, file), "n_global");
        REIS(1,
             fread(&entry_length_uint64, sizeof(entry_length_uint64), 1, file),
             "entry_length_uint64");
        nnode = (REF_GLOB)n_global;
        REIS(1, entry_length_uint64, "require entry_length == 1");

        REIS(1, fread(&count, sizeof(count), 1, file), "count");
        for (pair = 0; pair < count; pair++) {
          uint64_t length;
          /* key */
          REIS(1, fread(&length, sizeof(length), 1, file), "length");
          for (i = 0; i < (REF_INT)length; i++) {
            REIS(1, fread(&letter, sizeof(letter), 1, file), "number");
          }
          /* value */
          REIS(1, fread(&length, sizeof(length), 1, file), "length");
          for (i = 0; i < (REF_INT)length; i++) {
            REIS(1, fread(&letter, sizeof(letter), 1, file), "number");
          }
        }
      }

      if (verbose)
        printf("file nnode %ld ref nnode " REF_GLOB_FMT "\n", nnode,
               ref_node_n_global(ref_node));

      if ((nnode != ref_node_n_global(ref_node)) &&
          (nnode / 2 != ref_node_n_global(ref_node))) {
        printf("file " REF_GLOB_FMT " ref_node " REF_GLOB_FMT "\n", nnode,
               ref_node_n_global(ref_node));
        THROW("global count mismatch");
      }
    }
    RSS(ref_mpi_bcast(ref_node_mpi(ref_node), &nnode, 1, REF_INT_TYPE),
        "bcast ldim");

    chunk = (REF_INT)MAX(100000, nnode / ref_mpi_n(ref_node_mpi(ref_node)));
    chunk = (REF_INT)MIN((REF_GLOB)chunk, nnode);

    ref_malloc_init(data, chunk, REF_DBL, -1.0);

    nnode_read = 0;
    while (nnode_read < nnode) {
      if (verbose)
        printf("nnode " REF_GLOB_FMT " read " REF_GLOB_FMT "\n", nnode,
               nnode_read);

      section_size = (REF_INT)MIN((REF_GLOB)chunk, nnode - nnode_read);
      if (ref_mpi_once(ref_node_mpi(ref_node))) {
        REIS(section_size,
             fread(data, sizeof(REF_DBL), (size_t)(section_size), file), "dat");
        RSS(ref_mpi_bcast(ref_node_mpi(ref_node), data, chunk, REF_DBL_TYPE),
            "bcast");
      } else {
        RSS(ref_mpi_bcast(ref_node_mpi(ref_node), data, chunk, REF_DBL_TYPE),
            "bcast");
      }
      for (node = 0; node < section_size; node++) {
        global = node + nnode_read;
        RXS(ref_node_local(ref_node, global, &local), REF_NOT_FOUND, "local");
        if (REF_EMPTY != local) {
          (*scalar)[field + local * (*ldim)] = data[node];
        }
      }
      nnode_read += section_size;
    }

    ref_free(data);

    if (ref_mpi_once(ref_node_mpi(ref_node))) {
      if (nnode == ref_node_n_global(ref_node))
        REIS(next_position, ftello(file), "end location");
    }
  }

  if (ref_mpi_once(ref_node_mpi(ref_node))) {
    REIS(0, fclose(file), "close file");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_part_scalar(REF_GRID ref_grid, REF_INT *ldim, REF_DBL **scalar,
                           const char *filename) {
  REF_NODE ref_node = ref_grid_node(ref_grid);
  size_t end_of_string;

  end_of_string = strlen(filename);

  if (end_of_string > 12 &&
      strcmp(&filename[end_of_string - 12], ".restart_sol") == 0) {
    RSS(ref_part_scalar_cell_restart_sol(ref_grid, ldim, scalar, filename),
        "restart_sol failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 4 && strcmp(&filename[end_of_string - 4], ".rst") == 0) {
    RSS(ref_part_scalar_rst(ref_node, ldim, scalar, filename), "rst failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 4 && strcmp(&filename[end_of_string - 4], ".sol") == 0) {
    RSS(ref_part_scalar_sol(ref_node, ldim, scalar, filename), "sol failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 5 && strcmp(&filename[end_of_string - 5], ".solb") == 0) {
    RSS(ref_part_scalar_solb(ref_node, ldim, scalar, filename), "solb failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 5 && strcmp(&filename[end_of_string - 5], ".snap") == 0) {
    RSS(ref_part_scalar_snap(ref_node, ldim, scalar, filename), "snap failed");
    return REF_SUCCESS;
  }
  if (end_of_string > 4 && strcmp(&filename[end_of_string - 4], ".plt") == 0) {
    RSS(ref_part_scalar_plt(ref_grid, ldim, scalar, filename), "snap failed");
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
  if (end_of_string > 4 && strcmp(&filename[end_of_string - 4], ".avm") == 0) {
    RSS(ref_part_avm(ref_grid_ptr, ref_mpi, filename), "avm failed");
    return REF_SUCCESS;
  }
  printf("%s: %d: %s %s\n", __FILE__, __LINE__,
         "input file name extension unknown", filename);
  RSS(REF_FAILURE, "unknown file extension");
  return REF_FAILURE;
}
