
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

#include "ref_elast.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_grid.h"
#include "ref_malloc.h"
#include "ref_matrix.h"
#include "ref_node.h"

REF_STATUS ref_elast_create(REF_ELAST *ref_elast_ptr, REF_GRID ref_grid) {
  REF_ELAST ref_elast;

  ref_malloc(*ref_elast_ptr, 1, REF_ELAST_STRUCT);

  ref_elast = *ref_elast_ptr;

  ref_elast_grid(ref_elast) = ref_grid;
  RSS(ref_comprow_create(&(ref_elast->ref_comprow), ref_elast_grid(ref_elast)),
      "comprow");

  ref_malloc_init(ref_elast->a,
                  3 * 3 * ref_comprow_nnz(ref_elast_comprow(ref_elast)),
                  REF_DBL, 0.0);

  ref_malloc_init(ref_elast->displacement,
                  3 * ref_comprow_max(ref_elast_comprow(ref_elast)), REF_DBL,
                  0.0);
  ref_malloc_init(ref_elast->bc, ref_comprow_max(ref_elast_comprow(ref_elast)),
                  REF_INT, 0);

  return REF_SUCCESS;
}

REF_STATUS ref_elast_free(REF_ELAST ref_elast) {
  if (NULL == (void *)ref_elast) return REF_NULL;

  ref_free(ref_elast->bc);
  ref_free(ref_elast->displacement);
  ref_free(ref_elast->a);
  ref_comprow_free(ref_elast->ref_comprow);

  ref_free(ref_elast);

  return REF_SUCCESS;
}

REF_STATUS ref_elast_inspect(REF_ELAST ref_elast) {
  REF_COMPROW ref_comprow = ref_elast_comprow(ref_elast);
  REF_INT row, col, entry, i, j;
  for (row = 0; row < ref_comprow_max(ref_comprow); row++)
    if (ref_comprow->first[row + 1] > ref_comprow->first[row]) {
      for (entry = ref_comprow->first[row]; entry < ref_comprow->first[row + 1];
           entry++) {
        col = ref_comprow->col[entry];
        printf("row %d col %d\n", row, col);
        for (i = 0; i < 3; i++) {
          for (j = 0; j < 3; j++) {
            printf(" %10.5f", ref_elast->a[i + 3 * j + 9 * entry]);
          }
          printf("\n");
        }
      }
    }

  return REF_SUCCESS;
}

REF_STATUS ref_elast_displace(REF_ELAST ref_elast, REF_INT node,
                              REF_DBL *dxyz) {
  REF_GRID ref_grid = ref_elast_grid(ref_elast);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT i;
  if (!ref_node_valid(ref_node, node)) RSS(REF_INVALID, "invalid node index");
  for (i = 0; i < 3; i++) ref_elast->displacement[i + 3 * node] = dxyz[i];
  ref_elast->bc[node] = 1;
  return REF_SUCCESS;
}

REF_STATUS ref_elast_assemble(REF_ELAST ref_elast) {
  REF_COMPROW ref_comprow = ref_elast_comprow(ref_elast);
  REF_GRID ref_grid = ref_elast_grid(ref_elast);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT i;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node, entry;

  double x1, y1, z1;
  double x2, y2, z2;
  double x3, y3, z3;
  double x4, y4, z4;
  double nx1, ny1, nz1;
  double nx2, ny2, nz2;
  double nx3, ny3, nz3;
  double nx4, ny4, nz4;
  double vol;

  int idiag1, idiag2, idiag3, idiag4;
  int ioff12, ioff13, ioff14;
  int ioff21, ioff23, ioff24;
  int ioff31, ioff32, ioff34;
  int ioff41, ioff42, ioff43;

  double c;
  double c1x, c1y, c1z;
  double c2x, c2y, c2z;
  double c3x, c3y, c3z;
  double c4x, c4y, c4z;

  double xn, yn, zn;

  double aspect_ratio;
  double modulus_of_elasticity, poisson, mu, lambda;

  double xlambda, ylambda, zlambda;
  double xmu, ymu, zmu;
  double muterm;

  for (i = 0; i < 9 * ref_comprow_nnz(ref_comprow); i++) ref_elast->a[i] = 0.0;

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    x1 = ref_node_xyz(ref_node, 0, nodes[0]);
    y1 = ref_node_xyz(ref_node, 1, nodes[0]);
    z1 = ref_node_xyz(ref_node, 2, nodes[0]);
    x2 = ref_node_xyz(ref_node, 0, nodes[1]);
    y2 = ref_node_xyz(ref_node, 1, nodes[1]);
    z2 = ref_node_xyz(ref_node, 2, nodes[1]);
    x3 = ref_node_xyz(ref_node, 0, nodes[2]);
    y3 = ref_node_xyz(ref_node, 1, nodes[2]);
    z3 = ref_node_xyz(ref_node, 2, nodes[2]);
    x4 = ref_node_xyz(ref_node, 0, nodes[3]);
    y4 = ref_node_xyz(ref_node, 1, nodes[3]);
    z4 = ref_node_xyz(ref_node, 2, nodes[3]);

    nx1 = 0.5 * ((y2 - y4) * (z3 - z4) - (y3 - y4) * (z2 - z4));
    ny1 = 0.5 * ((z2 - z4) * (x3 - x4) - (z3 - z4) * (x2 - x4));
    nz1 = 0.5 * ((x2 - x4) * (y3 - y4) - (x3 - x4) * (y2 - y4));

    nx2 = 0.5 * ((y3 - y4) * (z1 - z4) - (y1 - y4) * (z3 - z4));
    ny2 = 0.5 * ((z3 - z4) * (x1 - x4) - (z1 - z4) * (x3 - x4));
    nz2 = 0.5 * ((x3 - x4) * (y1 - y4) - (x1 - x4) * (y3 - y4));

    nx3 = 0.5 * ((y1 - y4) * (z2 - z4) - (y2 - y4) * (z1 - z4));
    ny3 = 0.5 * ((z1 - z4) * (x2 - x4) - (z2 - z4) * (x1 - x4));
    nz3 = 0.5 * ((x1 - x4) * (y2 - y4) - (x2 - x4) * (y1 - y4));

    nx4 = -nx1 - nx2 - nx3;
    ny4 = -ny1 - ny2 - ny3;
    nz4 = -nz1 - nz2 - nz3;

    vol = (((y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1)) * (x4 - x1) -
           ((x2 - x1) * (z3 - z1) - (x3 - x1) * (z2 - z1)) * (y4 - y1) +
           ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) * (z4 - z1)) /
          6.0;

    if (vol <= 0.0) {
      printf("negative vol %f %f %f %e\n", x1, y1, z1, vol);
      RSS(REF_FAILURE, "neg vol");
    }

    RSS(ref_comprow_entry(ref_comprow, nodes[0], nodes[0], &idiag1), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[0], nodes[1], &ioff12), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[0], nodes[2], &ioff13), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[0], nodes[3], &ioff14), "e");

    RSS(ref_comprow_entry(ref_comprow, nodes[1], nodes[0], &ioff21), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[1], nodes[1], &idiag2), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[1], nodes[2], &ioff23), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[1], nodes[3], &ioff24), "e");

    RSS(ref_comprow_entry(ref_comprow, nodes[2], nodes[0], &ioff31), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[2], nodes[1], &ioff32), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[2], nodes[2], &idiag3), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[2], nodes[3], &ioff34), "e");

    RSS(ref_comprow_entry(ref_comprow, nodes[3], nodes[0], &ioff41), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[3], nodes[1], &ioff42), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[3], nodes[2], &ioff43), "e");
    RSS(ref_comprow_entry(ref_comprow, nodes[3], nodes[3], &idiag4), "e");

    idiag1 *= 9;
    ioff12 *= 9;
    ioff13 *= 9;
    ioff14 *= 9;
    ioff21 *= 9;
    idiag2 *= 9;
    ioff23 *= 9;
    ioff24 *= 9;
    ioff31 *= 9;
    ioff32 *= 9;
    idiag3 *= 9;
    ioff34 *= 9;
    ioff41 *= 9;
    ioff42 *= 9;
    ioff43 *= 9;
    idiag4 *= 9;

    c = 1.0 / (3.0 * vol);
    c1x = c * (nx2 + nx3 + nx4);
    c2x = c * (nx3 + nx4 + nx1);
    c3x = c * (nx4 + nx1 + nx2);
    c4x = c * (nx1 + nx2 + nx3);

    c1y = c * (ny2 + ny3 + ny4);
    c2y = c * (ny3 + ny4 + ny1);
    c3y = c * (ny4 + ny1 + ny2);
    c4y = c * (ny1 + ny2 + ny3);

    c1z = c * (nz2 + nz3 + nz4);
    c2z = c * (nz3 + nz4 + nz1);
    c3z = c * (nz4 + nz1 + nz2);
    c4z = c * (nz1 + nz2 + nz3);

    aspect_ratio = 1.0;
    modulus_of_elasticity = 1.0 / aspect_ratio;
    modulus_of_elasticity = pow(modulus_of_elasticity, 0.2);

    poisson = 0.0;

    mu = modulus_of_elasticity / 2.0 / (1.0 + poisson);
    lambda = poisson * modulus_of_elasticity / (1.0 + poisson) /
             (1.0 - 2.0 * poisson);

    /* node1 */
    xn = nx1 / 3.0;
    yn = ny1 / 3.0;
    zn = nz1 / 3.0;

    xlambda = xn * lambda;
    ylambda = yn * lambda;
    zlambda = zn * lambda;

    xmu = xn * mu;
    ymu = yn * mu;
    zmu = zn * mu;

    muterm = mu * (xn * c1x + yn * c1y + zn * c1z);

    ref_elast->a[0 + idiag1] += (muterm + xlambda * c1x + xmu * c1x);
    ref_elast->a[1 + idiag1] += (xlambda * c1y + ymu * c1x);
    ref_elast->a[2 + idiag1] += (xlambda * c1z + zmu * c1x);
    ref_elast->a[3 + idiag1] += (ylambda * c1x + xmu * c1y);
    ref_elast->a[4 + idiag1] += (muterm + ylambda * c1y + ymu * c1y);
    ref_elast->a[5 + idiag1] += (ylambda * c1z + zmu * c1y);
    ref_elast->a[6 + idiag1] += (zlambda * c1x + xmu * c1z);
    ref_elast->a[7 + idiag1] += (zlambda * c1y + ymu * c1z);
    ref_elast->a[8 + idiag1] += (muterm + zlambda * c1z + zmu * c1z);

    muterm = mu * (xn * c2x + yn * c2y + zn * c2z);

    ref_elast->a[0 + ioff12] += (muterm + xlambda * c2x + xmu * c2x);
    ref_elast->a[1 + ioff12] += (xlambda * c2y + ymu * c2x);
    ref_elast->a[2 + ioff12] += (xlambda * c2z + zmu * c2x);
    ref_elast->a[3 + ioff12] += (ylambda * c2x + xmu * c2y);
    ref_elast->a[4 + ioff12] += (muterm + ylambda * c2y + ymu * c2y);
    ref_elast->a[5 + ioff12] += (ylambda * c2z + zmu * c2y);
    ref_elast->a[6 + ioff12] += (zlambda * c2x + xmu * c2z);
    ref_elast->a[7 + ioff12] += (zlambda * c2y + ymu * c2z);
    ref_elast->a[8 + ioff12] += (muterm + zlambda * c2z + zmu * c2z);

    muterm = mu * (xn * c3x + yn * c3y + zn * c3z);

    ref_elast->a[0 + ioff13] += (muterm + xlambda * c3x + xmu * c3x);
    ref_elast->a[1 + ioff13] += (xlambda * c3y + ymu * c3x);
    ref_elast->a[2 + ioff13] += (xlambda * c3z + zmu * c3x);
    ref_elast->a[3 + ioff13] += (ylambda * c3x + xmu * c3y);
    ref_elast->a[4 + ioff13] += (muterm + ylambda * c3y + ymu * c3y);
    ref_elast->a[5 + ioff13] += (ylambda * c3z + zmu * c3y);
    ref_elast->a[6 + ioff13] += (zlambda * c3x + xmu * c3z);
    ref_elast->a[7 + ioff13] += (zlambda * c3y + ymu * c3z);
    ref_elast->a[8 + ioff13] += (muterm + zlambda * c3z + zmu * c3z);

    muterm = mu * (xn * c4x + yn * c4y + zn * c4z);

    ref_elast->a[0 + ioff14] += (muterm + xlambda * c4x + xmu * c4x);
    ref_elast->a[1 + ioff14] += (xlambda * c4y + ymu * c4x);
    ref_elast->a[2 + ioff14] += (xlambda * c4z + zmu * c4x);
    ref_elast->a[3 + ioff14] += (ylambda * c4x + xmu * c4y);
    ref_elast->a[4 + ioff14] += (muterm + ylambda * c4y + ymu * c4y);
    ref_elast->a[5 + ioff14] += (ylambda * c4z + zmu * c4y);
    ref_elast->a[6 + ioff14] += (zlambda * c4x + xmu * c4z);
    ref_elast->a[7 + ioff14] += (zlambda * c4y + ymu * c4z);
    ref_elast->a[8 + ioff14] += (muterm + zlambda * c4z + zmu * c4z);

    /* node2 */
    xn = nx2 / 3.0;
    yn = ny2 / 3.0;
    zn = nz2 / 3.0;

    xlambda = xn * lambda;
    ylambda = yn * lambda;
    zlambda = zn * lambda;

    xmu = xn * mu;
    ymu = yn * mu;
    zmu = zn * mu;

    muterm = mu * (xn * c2x + yn * c2y + zn * c2z);

    ref_elast->a[0 + idiag2] += (muterm + xlambda * c2x + xmu * c2x);
    ref_elast->a[1 + idiag2] += (xlambda * c2y + ymu * c2x);
    ref_elast->a[2 + idiag2] += (xlambda * c2z + zmu * c2x);
    ref_elast->a[3 + idiag2] += (ylambda * c2x + xmu * c2y);
    ref_elast->a[4 + idiag2] += (muterm + ylambda * c2y + ymu * c2y);
    ref_elast->a[5 + idiag2] += (ylambda * c2z + zmu * c2y);
    ref_elast->a[6 + idiag2] += (zlambda * c2x + xmu * c2z);
    ref_elast->a[7 + idiag2] += (zlambda * c2y + ymu * c2z);
    ref_elast->a[8 + idiag2] += (muterm + zlambda * c2z + zmu * c2z);

    muterm = mu * (xn * c3x + yn * c3y + zn * c3z);

    ref_elast->a[0 + ioff23] += (muterm + xlambda * c3x + xmu * c3x);
    ref_elast->a[1 + ioff23] += (xlambda * c3y + ymu * c3x);
    ref_elast->a[2 + ioff23] += (xlambda * c3z + zmu * c3x);
    ref_elast->a[3 + ioff23] += (ylambda * c3x + xmu * c3y);
    ref_elast->a[4 + ioff23] += (muterm + ylambda * c3y + ymu * c3y);
    ref_elast->a[5 + ioff23] += (ylambda * c3z + zmu * c3y);
    ref_elast->a[6 + ioff23] += (zlambda * c3x + xmu * c3z);
    ref_elast->a[7 + ioff23] += (zlambda * c3y + ymu * c3z);
    ref_elast->a[8 + ioff23] += (muterm + zlambda * c3z + zmu * c3z);

    muterm = mu * (xn * c4x + yn * c4y + zn * c4z);

    ref_elast->a[0 + ioff24] += (muterm + xlambda * c4x + xmu * c4x);
    ref_elast->a[1 + ioff24] += (xlambda * c4y + ymu * c4x);
    ref_elast->a[2 + ioff24] += (xlambda * c4z + zmu * c4x);
    ref_elast->a[3 + ioff24] += (ylambda * c4x + xmu * c4y);
    ref_elast->a[4 + ioff24] += (muterm + ylambda * c4y + ymu * c4y);
    ref_elast->a[5 + ioff24] += (ylambda * c4z + zmu * c4y);
    ref_elast->a[6 + ioff24] += (zlambda * c4x + xmu * c4z);
    ref_elast->a[7 + ioff24] += (zlambda * c4y + ymu * c4z);
    ref_elast->a[8 + ioff24] += (muterm + zlambda * c4z + zmu * c4z);

    muterm = mu * (xn * c1x + yn * c1y + zn * c1z);

    ref_elast->a[0 + ioff21] += (muterm + xlambda * c1x + xmu * c1x);
    ref_elast->a[1 + ioff21] += (xlambda * c1y + ymu * c1x);
    ref_elast->a[2 + ioff21] += (xlambda * c1z + zmu * c1x);
    ref_elast->a[3 + ioff21] += (ylambda * c1x + xmu * c1y);
    ref_elast->a[4 + ioff21] += (muterm + ylambda * c1y + ymu * c1y);
    ref_elast->a[5 + ioff21] += (ylambda * c1z + zmu * c1y);
    ref_elast->a[6 + ioff21] += (zlambda * c1x + xmu * c1z);
    ref_elast->a[7 + ioff21] += (zlambda * c1y + ymu * c1z);
    ref_elast->a[8 + ioff21] += (muterm + zlambda * c1z + zmu * c1z);

    /* node3 */
    xn = nx3 / 3.0;
    yn = ny3 / 3.0;
    zn = nz3 / 3.0;

    xlambda = xn * lambda;
    ylambda = yn * lambda;
    zlambda = zn * lambda;

    xmu = xn * mu;
    ymu = yn * mu;
    zmu = zn * mu;

    muterm = mu * (xn * c3x + yn * c3y + zn * c3z);

    ref_elast->a[0 + idiag3] += (muterm + xlambda * c3x + xmu * c3x);
    ref_elast->a[1 + idiag3] += (xlambda * c3y + ymu * c3x);
    ref_elast->a[2 + idiag3] += (xlambda * c3z + zmu * c3x);
    ref_elast->a[3 + idiag3] += (ylambda * c3x + xmu * c3y);
    ref_elast->a[4 + idiag3] += (muterm + ylambda * c3y + ymu * c3y);
    ref_elast->a[5 + idiag3] += (ylambda * c3z + zmu * c3y);
    ref_elast->a[6 + idiag3] += (zlambda * c3x + xmu * c3z);
    ref_elast->a[7 + idiag3] += (zlambda * c3y + ymu * c3z);
    ref_elast->a[8 + idiag3] += (muterm + zlambda * c3z + zmu * c3z);

    muterm = mu * (xn * c4x + yn * c4y + zn * c4z);

    ref_elast->a[0 + ioff34] += (muterm + xlambda * c4x + xmu * c4x);
    ref_elast->a[1 + ioff34] += (xlambda * c4y + ymu * c4x);
    ref_elast->a[2 + ioff34] += (xlambda * c4z + zmu * c4x);
    ref_elast->a[3 + ioff34] += (ylambda * c4x + xmu * c4y);
    ref_elast->a[4 + ioff34] += (muterm + ylambda * c4y + ymu * c4y);
    ref_elast->a[5 + ioff34] += (ylambda * c4z + zmu * c4y);
    ref_elast->a[6 + ioff34] += (zlambda * c4x + xmu * c4z);
    ref_elast->a[7 + ioff34] += (zlambda * c4y + ymu * c4z);
    ref_elast->a[8 + ioff34] += (muterm + zlambda * c4z + zmu * c4z);

    muterm = mu * (xn * c1x + yn * c1y + zn * c1z);

    ref_elast->a[0 + ioff31] += (muterm + xlambda * c1x + xmu * c1x);
    ref_elast->a[1 + ioff31] += (xlambda * c1y + ymu * c1x);
    ref_elast->a[2 + ioff31] += (xlambda * c1z + zmu * c1x);
    ref_elast->a[3 + ioff31] += (ylambda * c1x + xmu * c1y);
    ref_elast->a[4 + ioff31] += (muterm + ylambda * c1y + ymu * c1y);
    ref_elast->a[5 + ioff31] += (ylambda * c1z + zmu * c1y);
    ref_elast->a[6 + ioff31] += (zlambda * c1x + xmu * c1z);
    ref_elast->a[7 + ioff31] += (zlambda * c1y + ymu * c1z);
    ref_elast->a[8 + ioff31] += (muterm + zlambda * c1z + zmu * c1z);

    muterm = mu * (xn * c2x + yn * c2y + zn * c2z);

    ref_elast->a[0 + ioff32] += (muterm + xlambda * c2x + xmu * c2x);
    ref_elast->a[1 + ioff32] += (xlambda * c2y + ymu * c2x);
    ref_elast->a[2 + ioff32] += (xlambda * c2z + zmu * c2x);
    ref_elast->a[3 + ioff32] += (ylambda * c2x + xmu * c2y);
    ref_elast->a[4 + ioff32] += (muterm + ylambda * c2y + ymu * c2y);
    ref_elast->a[5 + ioff32] += (ylambda * c2z + zmu * c2y);
    ref_elast->a[6 + ioff32] += (zlambda * c2x + xmu * c2z);
    ref_elast->a[7 + ioff32] += (zlambda * c2y + ymu * c2z);
    ref_elast->a[8 + ioff32] += (muterm + zlambda * c2z + zmu * c2z);

    /* node4 */
    xn = nx4 / 3.0;
    yn = ny4 / 3.0;
    zn = nz4 / 3.0;

    xlambda = xn * lambda;
    ylambda = yn * lambda;
    zlambda = zn * lambda;

    xmu = xn * mu;
    ymu = yn * mu;
    zmu = zn * mu;

    muterm = mu * (xn * c4x + yn * c4y + zn * c4z);

    ref_elast->a[0 + idiag4] += (muterm + xlambda * c4x + xmu * c4x);
    ref_elast->a[1 + idiag4] += (xlambda * c4y + ymu * c4x);
    ref_elast->a[2 + idiag4] += (xlambda * c4z + zmu * c4x);
    ref_elast->a[3 + idiag4] += (ylambda * c4x + xmu * c4y);
    ref_elast->a[4 + idiag4] += (muterm + ylambda * c4y + ymu * c4y);
    ref_elast->a[5 + idiag4] += (ylambda * c4z + zmu * c4y);
    ref_elast->a[6 + idiag4] += (zlambda * c4x + xmu * c4z);
    ref_elast->a[7 + idiag4] += (zlambda * c4y + ymu * c4z);
    ref_elast->a[8 + idiag4] += (muterm + zlambda * c4z + zmu * c4z);

    muterm = mu * (xn * c1x + yn * c1y + zn * c1z);

    ref_elast->a[0 + ioff41] += (muterm + xlambda * c1x + xmu * c1x);
    ref_elast->a[1 + ioff41] += (xlambda * c1y + ymu * c1x);
    ref_elast->a[2 + ioff41] += (xlambda * c1z + zmu * c1x);
    ref_elast->a[3 + ioff41] += (ylambda * c1x + xmu * c1y);
    ref_elast->a[4 + ioff41] += (muterm + ylambda * c1y + ymu * c1y);
    ref_elast->a[5 + ioff41] += (ylambda * c1z + zmu * c1y);
    ref_elast->a[6 + ioff41] += (zlambda * c1x + xmu * c1z);
    ref_elast->a[7 + ioff41] += (zlambda * c1y + ymu * c1z);
    ref_elast->a[8 + ioff41] += (muterm + zlambda * c1z + zmu * c1z);

    muterm = mu * (xn * c2x + yn * c2y + zn * c2z);

    ref_elast->a[0 + ioff42] += (muterm + xlambda * c2x + xmu * c2x);
    ref_elast->a[1 + ioff42] += (xlambda * c2y + ymu * c2x);
    ref_elast->a[2 + ioff42] += (xlambda * c2z + zmu * c2x);
    ref_elast->a[3 + ioff42] += (ylambda * c2x + xmu * c2y);
    ref_elast->a[4 + ioff42] += (muterm + ylambda * c2y + ymu * c2y);
    ref_elast->a[5 + ioff42] += (ylambda * c2z + zmu * c2y);
    ref_elast->a[6 + ioff42] += (zlambda * c2x + xmu * c2z);
    ref_elast->a[7 + ioff42] += (zlambda * c2y + ymu * c2z);
    ref_elast->a[8 + ioff42] += (muterm + zlambda * c2z + zmu * c2z);

    muterm = mu * (xn * c3x + yn * c3y + zn * c3z);

    ref_elast->a[0 + ioff43] += (muterm + xlambda * c3x + xmu * c3x);
    ref_elast->a[1 + ioff43] += (xlambda * c3y + ymu * c3x);
    ref_elast->a[2 + ioff43] += (xlambda * c3z + zmu * c3x);
    ref_elast->a[3 + ioff43] += (ylambda * c3x + xmu * c3y);
    ref_elast->a[4 + ioff43] += (muterm + ylambda * c3y + ymu * c3y);
    ref_elast->a[5 + ioff43] += (ylambda * c3z + zmu * c3y);
    ref_elast->a[6 + ioff43] += (zlambda * c3x + xmu * c3z);
    ref_elast->a[7 + ioff43] += (zlambda * c3y + ymu * c3z);
    ref_elast->a[8 + ioff43] += (muterm + zlambda * c3z + zmu * c3z);
  }

  each_ref_node_valid_node(ref_node, node) {
    if (!ref_node_owned(ref_node, node) || 1 == ref_elast->bc[node]) {
      each_ref_comprow_row_entry(ref_comprow, node, entry) {
        for (i = 0; i < 9; i++) ref_elast->a[i + 9 * entry] = 0.0;
      }
      RSS(ref_comprow_entry(ref_comprow, node, node, &entry), "e");
      ref_elast->a[0 + 9 * entry] = 1.0;
      ref_elast->a[4 + 9 * entry] = 1.0;
      ref_elast->a[8 + 9 * entry] = 1.0;
    }
  }

  /* to set ghost node bcs */
  RSS(ref_node_ghost_int(ref_node, ref_elast->bc, 1), "ghost bcs");
  RSS(ref_node_ghost_dbl(ref_node, ref_elast->displacement, 3), "ghost disp");

  return REF_SUCCESS;
}

REF_STATUS ref_elast_relax(REF_ELAST ref_elast, REF_DBL *l2norm) {
  REF_COMPROW ref_comprow = ref_elast_comprow(ref_elast);
  REF_GRID ref_grid = ref_elast_grid(ref_elast);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_MPI ref_mpi = ref_grid_mpi(ref_grid);
  double ab[12];
  int entry, row, col, i, j;

  *l2norm = 0.0;
  each_ref_node_valid_node(ref_node, row) {
    if (ref_node_owned(ref_node, row) && 0 == ref_elast->bc[row]) {
      ab[9] = ab[10] = ab[11] = 0.0;
      each_ref_comprow_row_entry(ref_comprow, row, entry) {
        col = ref_comprow->col[entry];
        if (row != col) {
          for (i = 0; i < 3; i++) {
            ab[9 + i] -= (ref_elast->a[i + 0 * 3 + 9 * entry] *
                              ref_elast->displacement[0 + 3 * col] +
                          ref_elast->a[i + 1 * 3 + 9 * entry] *
                              ref_elast->displacement[1 + 3 * col] +
                          ref_elast->a[i + 2 * 3 + 9 * entry] *
                              ref_elast->displacement[2 + 3 * col]);
          }
        }
      }
      RSS(ref_comprow_entry(ref_comprow, row, row, &entry), "diag");
      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
          ab[i + 3 * j] = ref_elast->a[i + j * 3 + 9 * entry];
      RSS(ref_matrix_solve_ab(3, 4, ab), "solve");
      for (i = 0; i < 3; i++)
        *l2norm += pow(ref_elast->displacement[i + 3 * row] - ab[9 + i], 2);
      for (i = 0; i < 3; i++) ref_elast->displacement[i + 3 * row] = ab[9 + i];
    }
  }
  RSS(ref_node_ghost_dbl(ref_node, ref_elast->displacement, 3), "ghost disp");
  RSS(ref_mpi_allsum(ref_mpi, l2norm, 1, REF_DBL_TYPE),
      "sum l2 norm over parts");
  *l2norm /= (REF_DBL)ref_node_n_global(ref_node);
  *l2norm = sqrt(*l2norm);

  return REF_SUCCESS;
}
