
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

#include "ref_inflate.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_adj.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_sort.h"

int main() {
  {
    REF_DBL H, h0, rate;
    REF_INT n;
    H = 10.0;
    h0 = 1.0;
    n = 10;
    RSS(ref_inflate_rate(n, h0, H, &rate), "rate");
    RWDS(1.0, rate, 1.0e-8, "uniform");
  }

  {
    REF_DBL H, h0, rate;
    REF_INT n;
    H = 10.0;
    h0 = 0.5;
    n = 10;
    RSS(ref_inflate_rate(n, h0, H, &rate), "rate");
    RWDS(1.1469, rate, 1.0e-3, "growth");
  }

  /*
  SKIP_BLOCK("intel compiler") {
    REF_DBL H, h0, rate;
    REF_INT n;
    H = 0.1;
    h0 = 0.01;
    n = 10;
    RSS(ref_inflate_rate(n, h0, H, &rate), "rate");
    RWDS(1.0, rate, 1.0e-6, "growth");
  }
  */

  {
    REF_DBL H, h0, rate;
    REF_INT n;
    h0 = 0.5;
    n = 10;
    rate = 1.1469;
    RSS(ref_inflate_total_thickness(n, h0, rate, &H), "total");
    RWDS(10.0, H, 1.0e-3, "not total");
  }

  {
    REF_DBL dHdr, h0, rate;
    REF_INT n;
    REF_DBL tol = 1.0e-6, fd, Hp, Hm;
    h0 = 0.5;
    n = 10;
    rate = 1.1469;
    RSS(ref_inflate_total_thickness(n, h0, rate + tol, &Hp), "total+");
    RSS(ref_inflate_total_thickness(n, h0, rate - tol, &Hm), "total-");
    fd = (Hp - Hm) / 2.0 / tol;

    RSS(ref_inflate_dthickness(n, h0, rate, &dHdr), "deriv");

    RWDS(fd, dHdr, 1.0e-8, "not total");
  }

  {
    char file[] = "ref_inflate_test.mapbc";
    char family[] = "cigar";
    REF_INT bc = 3;
    FILE *f;
    REF_DICT ref_dict;
    REF_INT id, type;

    f = fopen(file, "w");
    fprintf(f, "#\n#\n#\n#\n");
    fprintf(f, "1219            3       3        0 0 SURFACE.  \n");
    fprintf(f, "1220            4       4        0 0 SURFACE.  \n");
    fprintf(f, "1221            2       2        0 0 cigar      \n");
    fprintf(f, "1222            0       0        0 0 cigar     \n");
    fprintf(f, "1223            3       3        0 0 cigar\n");
    fprintf(f, "1224            3       3        0 0 cigar       \n");
    fclose(f);

    RSS(ref_dict_create(&ref_dict), "create");

    RSS(ref_inflate_read_usm3d_mapbc(ref_dict, file, family, bc), "read mapbc");

    REIS(2, ref_dict_n(ref_dict), "lines");
    id = 1223;
    RSS(ref_dict_value(ref_dict, id, &type), "retrieve");
    REIS(REF_EMPTY, type, "type");
    id = 1224;
    RSS(ref_dict_value(ref_dict, id, &type), "retrieve");
    REIS(REF_EMPTY, type, "type");

    RSS(ref_dict_free(ref_dict), "free");
    REIS(0, remove(file), "test clean up");
  }

  return 0;
}
