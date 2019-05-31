
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
#include <string.h>

#include "ref_fortran.h"

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_export.h"
#include "ref_grid.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_node.h"

#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_mpi.h"
#include "ref_sort.h"
#include "ref_subdiv.h"

#include "ref_adapt.h"
#include "ref_metric.h"
#include "ref_migrate.h"
#include "ref_validation.h"

#include "ref_collapse.h"
#include "ref_face.h"
#include "ref_smooth.h"
#include "ref_split.h"
#include "ref_twod.h"

#include "ref_gather.h"
#include "ref_histogram.h"

int main(int argc, char *argv[]) {
  REF_INT nnodes, nnodes0;
  REF_INT nnodesg;
  REF_INT *l2g;
  REF_INT *part;
  REF_INT partition;
  REF_DBL *x;
  REF_DBL *y;
  REF_DBL *z;
  REF_INT *c2n;
  REF_INT *f2n;
  REF_DBL *m;
  REF_DBL *ratio;

  REF_DBL *aux;
  REF_INT naux, offset;

  REF_INT node_per_cell;
  REF_INT ncell;
  REF_INT node_per_face;
  REF_INT nface;
  REF_INT ibound;

  REF_INT node;

  RSS(ref_mpi_start(argc, argv), "start");

  nnodes = 4;
  nnodesg = 4;

  l2g = (REF_INT *)malloc(sizeof(REF_INT) * (size_t)nnodes);
  part = (REF_INT *)malloc(sizeof(REF_INT) * (size_t)nnodes);
  x = (REF_DBL *)malloc(sizeof(REF_DBL) * (size_t)nnodes);
  y = (REF_DBL *)malloc(sizeof(REF_DBL) * (size_t)nnodes);
  z = (REF_DBL *)malloc(sizeof(REF_DBL) * (size_t)nnodes);
  m = (REF_DBL *)malloc(sizeof(REF_DBL) * 6 * (size_t)nnodes);
  ratio = (REF_DBL *)malloc(sizeof(REF_DBL) * (size_t)nnodes);

  naux = 2;
  aux = (REF_DBL *)malloc(sizeof(REF_DBL) * (size_t)naux * (size_t)nnodes);

  l2g[0] = 3;
  part[0] = 1;
  x[0] = 0;
  y[0] = 1;
  z[0] = 0;
  l2g[1] = 4;
  part[1] = 1;
  x[1] = 0;
  y[1] = 0;
  z[1] = 1;
  l2g[2] = 1;
  part[2] = 2;
  x[2] = 0;
  y[2] = 0;
  z[2] = 0;
  l2g[3] = 2;
  part[3] = 2;
  x[3] = 1;
  y[3] = 0;
  z[3] = 0;

  for (node = 0; node < nnodes; node++) {
    m[0 + 6 * node] = 1.0;
    m[1 + 6 * node] = 0.0;
    m[2 + 6 * node] = 0.0;
    m[3 + 6 * node] = 1.0;
    m[4 + 6 * node] = 0.0;
    m[5 + 6 * node] = 1.0;
    ratio[node] = 1.0;
    aux[0 + naux * node] = 2.0;
    aux[1 + naux * node] = 4.0;
  }

  partition = 0;

  ref_fortran_allow_screen_output = REF_FALSE;

  RSS(REF_FORT_(ref_fortran_init, REF_FORTRAN_INIT)(&nnodes, &nnodesg, l2g,
                                                    part, &partition, x, y, z),
      "init node");

  node_per_cell = 4;
  ncell = 1;
  c2n = (REF_INT *)malloc(sizeof(REF_INT) * (size_t)node_per_cell *
                          (size_t)ncell);

  c2n[0] = 1;
  c2n[1] = 2;
  c2n[2] = 3;
  c2n[3] = 4;

  RSS(REF_FORT_(ref_fortran_import_cell, REF_FORTRAN_IMPORT_CELL)(
          &node_per_cell, &ncell, c2n),
      "import cell");

  node_per_face = 3;
  nface = 1;
  f2n = (REF_INT *)malloc(sizeof(REF_INT) * (size_t)node_per_face *
                          (size_t)nface);
  f2n[0] = 3;
  f2n[1] = 4;
  f2n[2] = 1;

  ibound = 1;
  RSS(REF_FORT_(ref_fortran_import_face, REF_FORTRAN_IMPORT_FACE)(
          &ibound, &node_per_face, &nface, f2n),
      "import face");

  RSS(REF_FORT_(ref_fortran_import_metric, REF_FORTRAN_IMPORT_METRIC)(&nnodes,
                                                                      m),
      "import metric");

  RSS(REF_FORT_(ref_fortran_import_ratio, REF_FORTRAN_IMPORT_RATIO)(&nnodes,
                                                                    ratio),
      "import ratio");

  RSS(REF_FORT_(ref_fortran_naux, REF_FORTRAN_NAUX)(&naux), "naux");

  offset = 0;
  RSS(REF_FORT_(ref_fortran_import_aux, REF_FORTRAN_IMPORT_AUX)(&naux, &nnodes,
                                                                &offset, aux),
      "aux");

  free(f2n);
  free(c2n);
  free(aux);
  free(ratio);
  free(m);
  free(z);
  free(y);
  free(x);
  free(part);
  free(l2g);

  nnodes = nnodesg = REF_EMPTY;
  RSS(REF_FORT_(ref_fortran_size_node, REF_FORTRAN_SIZE_NODE)(&nnodes0, &nnodes,
                                                              &nnodesg),
      "size_node");
  REIS(2, nnodes0, "n");
  REIS(4, nnodes, "n");
  REIS(4, nnodesg, "n");

  l2g = (REF_INT *)malloc(sizeof(REF_INT) * (size_t)nnodes);
  x = (REF_DBL *)malloc(sizeof(REF_DBL) * (size_t)nnodes);
  y = (REF_DBL *)malloc(sizeof(REF_DBL) * (size_t)nnodes);
  z = (REF_DBL *)malloc(sizeof(REF_DBL) * (size_t)nnodes);

  RSS(REF_FORT_(ref_fortran_node, REF_FORTRAN_NODE)(&nnodes, l2g, x, y, z),
      "get node");

  aux = (REF_DBL *)malloc(sizeof(REF_DBL) * (size_t)naux * (size_t)nnodes);

  RSS(REF_FORT_(ref_fortran_aux, REF_FORTRAN_NODE)(&naux, &nnodes, &offset,
                                                   aux),
      "get aux");

  ncell = REF_EMPTY;
  RSS(REF_FORT_(ref_fortran_size_cell, REF_FORTRAN_SIZE_CELL)(&node_per_cell,
                                                              &ncell),
      "size cell");
  REIS(1, ncell, "n");

  c2n = (REF_INT *)malloc(sizeof(REF_INT) * (size_t)node_per_cell *
                          (size_t)ncell);

  RSS(REF_FORT_(ref_fortran_cell, REF_FORTRAN_CELL)(&node_per_cell, &ncell,
                                                    c2n),
      "get cell");

  ibound = 1;
  node_per_face = 3;
  nface = REF_EMPTY;
  RSS(REF_FORT_(ref_fortran_size_face, REF_FORTRAN_SIZE_FACE)(
          &ibound, &node_per_face, &nface),
      "size face");
  REIS(1, nface, "n");

  f2n = (REF_INT *)malloc(sizeof(REF_INT) * (size_t)node_per_face *
                          (size_t)nface);

  RSS(REF_FORT_(ref_fortran_face, REF_FORTRAN_FACE)(&ibound, &node_per_face,
                                                    &nface, f2n),
      "face");

  free(f2n);
  free(c2n);
  free(aux);
  free(z);
  free(y);
  free(x);
  free(l2g);

  RSS(REF_FORT_(ref_fortran_free, REF_FORTRAN_FREE)(), "free");

  RSS(ref_mpi_stop(), "stop");

  return 0;
}
