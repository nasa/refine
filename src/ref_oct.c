
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

#include "ref_malloc.h"

REF_FCN REF_STATUS ref_oct_create(REF_OCT *ref_oct_ptr) {
  REF_OCT ref_oct;

  ref_malloc(*ref_oct_ptr, 1, REF_OCT_STRUCT);
  ref_oct = (*ref_oct_ptr);

  ref_oct->bbox[0] = 0.0;
  ref_oct->bbox[1] = 1.0;
  ref_oct->bbox[2] = 0.0;
  ref_oct->bbox[3] = 1.0;
  ref_oct->bbox[4] = 0.0;
  ref_oct->bbox[5] = 1.0;

  ref_oct->n = 1;
  ref_oct->max = 9;
  ref_oct->children = NULL;
  ref_malloc_init(ref_oct->children, 8 * ref_oct->max, REF_INT, REF_EMPTY);

  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_free(REF_OCT ref_oct) {
  if (NULL == (void *)ref_oct) return REF_NULL;
  ref_free(ref_oct->children);
  ref_free(ref_oct);
  return REF_SUCCESS;
}

REF_FCN REF_STATUS ref_oct_tec(REF_OCT ref_oct, const char *filename) {
  FILE *f;
  const char *zonetype = "febrick";
  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  fprintf(f, "title=\"tecplot refine octree\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\"\n");

  fprintf(f,
          "zone t=\"octree\", nodes=%d, elements=%d, datapacking=%s, "
          "zonetype=%s\n",
          8, 1, "point", zonetype);

  fprintf(f, "%f %f %f\n", ref_oct->bbox[0], ref_oct->bbox[2],
          ref_oct->bbox[4]);
  fprintf(f, "%f %f %f\n", ref_oct->bbox[1], ref_oct->bbox[2],
          ref_oct->bbox[4]);
  fprintf(f, "%f %f %f\n", ref_oct->bbox[1], ref_oct->bbox[3],
          ref_oct->bbox[4]);
  fprintf(f, "%f %f %f\n", ref_oct->bbox[0], ref_oct->bbox[3],
          ref_oct->bbox[4]);
  fprintf(f, "%f %f %f\n", ref_oct->bbox[0], ref_oct->bbox[2],
          ref_oct->bbox[5]);
  fprintf(f, "%f %f %f\n", ref_oct->bbox[1], ref_oct->bbox[2],
          ref_oct->bbox[5]);
  fprintf(f, "%f %f %f\n", ref_oct->bbox[1], ref_oct->bbox[3],
          ref_oct->bbox[5]);
  fprintf(f, "%f %f %f\n", ref_oct->bbox[0], ref_oct->bbox[3],
          ref_oct->bbox[5]);

  fprintf(f, "%d %d %d %d %d %d %d %d\n", 1, 2, 3, 4, 5, 6, 7, 8);

  fclose(f);
  return REF_SUCCESS;
}
