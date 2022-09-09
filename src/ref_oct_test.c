
/* Copyright 2006, 2014, 2021 United States Government as represented
 * by the Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine version 3 unstructured grid adaptation platform is
 * licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * https://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#include "ref_oct.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_args.h"
#include "ref_mpi.h"

int main(int argc, char *argv[]) {
  REF_INT pos;
  REF_MPI ref_mpi;
  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "make mpi");

  RXS(ref_args_find(argc, argv, "--tec", &pos), REF_NOT_FOUND, "arg search");
  if (pos != REF_EMPTY && pos + 1 < argc) {
    REF_OCT ref_oct;
    RSS(ref_oct_create(&ref_oct), "make oct");
    RSS(ref_oct_split(ref_oct, 0), "split root");
    RSS(ref_oct_split(ref_oct, 1), "split first child");
    RSS(ref_oct_split(ref_oct, 9), "split second gen");
    RSS(ref_oct_tec(ref_oct, argv[pos + 1]), "tec");
    RSS(ref_oct_free(ref_oct), "search oct");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  RXS(ref_args_find(argc, argv, "--to", &pos), REF_NOT_FOUND, "arg search");
  if (pos != REF_EMPTY && pos + 1 < argc) {
    REF_OCT ref_oct;
    REF_DBL xyz[3], h;
    RSS(ref_oct_create(&ref_oct), "make oct");
    xyz[0] = 0.3;
    xyz[1] = 0.4;
    xyz[2] = 0.5;
    h = 0.01;
    RSS(ref_oct_split_at(ref_oct, xyz, h), "split xyz h");
    RSS(ref_oct_tec(ref_oct, argv[pos + 1]), "tec");
    RSS(ref_oct_free(ref_oct), "search oct");
    RSS(ref_mpi_free(ref_mpi), "mpi free");
    RSS(ref_mpi_stop(), "stop");
    return 0;
  }

  { /* create */
    REF_OCT ref_oct;
    RSS(ref_oct_create(&ref_oct), "make oct");
    RSS(ref_oct_free(ref_oct), "free oct");
  }

  { /* split root */
    REF_OCT ref_oct;
    RSS(ref_oct_create(&ref_oct), "make oct");
    RSS(ref_oct_split(ref_oct, 0), "split oct");
    RSS(ref_oct_free(ref_oct), "free oct");
  }

  { /* child's bbox */
    REF_DBL parent_bbox[6], child_bbox[6];
    REF_INT child_index = 0;
    REF_DBL tol = -1.0;
    parent_bbox[0] = 0.0;
    parent_bbox[1] = 1.0;
    parent_bbox[2] = 0.0;
    parent_bbox[3] = 1.0;
    parent_bbox[4] = 0.0;
    parent_bbox[5] = 1.0;
    RSS(ref_oct_child_bbox(parent_bbox, child_index, child_bbox), "bbox");
    RWDS(0.0, child_bbox[0], tol, "not zero");
  }

  { /* bbox diag */
    REF_DBL bbox[6];
    REF_DBL diag, tol = -1.0;
    bbox[0] = 0.0;
    bbox[1] = 0.5;
    bbox[2] = 0.5;
    bbox[3] = 1.0;
    bbox[4] = 0.2;
    bbox[5] = 0.7;
    RSS(ref_oct_bbox_diag(bbox, &diag), "bbox diag");
    RWDS(sqrt(0.75), diag, tol, "wrong size");
  }

  { /* contains root */
    REF_OCT ref_oct;
    REF_DBL xyz[] = {0.1, 0.1, 0.1};
    REF_DBL bbox[6];
    REF_INT node;
    RSS(ref_oct_create(&ref_oct), "make oct");
    RSS(ref_oct_contains(ref_oct, xyz, &node, bbox), "contains oct");
    REIS(0, node, "expects root");
    xyz[0] = 100.0;
    RSS(ref_oct_contains(ref_oct, xyz, &node, bbox), "contains oct");
    REIS(REF_EMPTY, node, "expects empty outside root");
    RSS(ref_oct_free(ref_oct), "free oct");
  }

  { /* contains child */
    REF_OCT ref_oct;
    REF_DBL xyz[] = {0.1, 0.1, 0.1};
    REF_DBL bbox[6];
    REF_INT node;
    RSS(ref_oct_create(&ref_oct), "make oct");
    RSS(ref_oct_split(ref_oct, 0), "split oct");
    RSS(ref_oct_contains(ref_oct, xyz, &node, bbox), "contains oct");
    REIS(1, node, "expects first child");
    xyz[0] = 0.9;
    xyz[1] = 0.1;
    xyz[2] = 0.1;
    RSS(ref_oct_contains(ref_oct, xyz, &node, bbox), "contains oct");
    REIS(2, node, "expects last child");
    xyz[0] = 0.1;
    xyz[1] = 0.9;
    xyz[2] = 0.1;
    RSS(ref_oct_contains(ref_oct, xyz, &node, bbox), "contains oct");
    REIS(4, node, "expects last child");
    xyz[0] = 0.9;
    xyz[1] = 0.9;
    xyz[2] = 0.9;
    RSS(ref_oct_contains(ref_oct, xyz, &node, bbox), "contains oct");
    REIS(7, node, "expects last child");
    RSS(ref_oct_free(ref_oct), "free oct");
  }

  RSS(ref_mpi_free(ref_mpi), "mpi free");
  RSS(ref_mpi_stop(), "stop");
  return 0;
}
