
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

#include "ref_html.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ref_malloc.h"
#include "ref_math.h"
#include "ref_matrix.h"

REF_STATUS ref_html_create(REF_HTML *ref_html_ptr, const char *filename) {
  FILE *f;

  f = fopen(filename, "w");
  if (NULL == (void *)f) printf("unable to open %s\n", filename);
  RNS(f, "unable to open file");

  ref_malloc(*ref_html_ptr, 1, REF_HTML_STRUCT);
  ref_html_file(*ref_html_ptr) = f;

  fprintf(f, "<html>\n");
  fprintf(f, "  <head>\n");
  fprintf(f, "    <title>refine export</title>\n");
  fprintf(f, "    <link rel='stylesheet' type='text/css'\n");
  fprintf(f, "          href='http://www.x3dom.org/download/x3dom.css'>\n");
  fprintf(f, "    </link>\n");
  fprintf(f, "    <script type='text/javascript'\n");
  fprintf(f, "          src='http://www.x3dom.org/download/x3dom.js'>\n");
  fprintf(f, "    </script>\n");
  fprintf(f, "    <style>\n");
  fprintf(f, "      x3d {width:100%%;height:100%%;border:none}\n");
  fprintf(f, "      body {margin:0;width:100%%;height:100%%;}\n");
  fprintf(f, "    </style>\n");
  fprintf(f, "  </head>\n");
  fprintf(f, "  <body id='body'>\n");
  fprintf(f, "    <a href=\"http://x3dom.org/docs/dev/navigation.html\">\n");
  fprintf(f, "camera control help\n");
  fprintf(f, "    </a>\n");
  fprintf(f, "    <x3d id='x3d'><scene><shape>\n");

  return REF_SUCCESS;
}

/* http://www.web3d.org/x3d/content/examples/Conformance/Geometry/IndexedLineSet/_pages/page09.html
 */

REF_STATUS ref_html_diagonal_system(REF_HTML ref_html, REF_DBL *origin,
                                    REF_DBL *d) {
  FILE *f = ref_html_file(ref_html);
  REF_INT i, j, n;
  REF_DBL dt, t;
  REF_DBL x, y, z;
  n = 72;
  dt = 2 * ref_math_pi / (REF_DBL)n;

  fprintf(f, "      <IndexedLineSet coordIndex='\n");
  for (j = 0; j < 3; j++) {
    for (i = 0; i < n; i++) fprintf(f, " %d", j * n + i);
    fprintf(f, " %d %d\n", j * n + 0, -1);
  }
  fprintf(f, "      ' >\n");

  fprintf(f, "      <Coordinate point='\n");
  for (i = 0; i < n; i++) {
    t = dt * (REF_DBL)i;

    x = origin[0] + ref_matrix_eig(d, 0) * ref_matrix_vec(d, 0, 0) * sin(t) +
        ref_matrix_eig(d, 1) * ref_matrix_vec(d, 0, 1) * cos(t);

    y = origin[1] + ref_matrix_eig(d, 0) * ref_matrix_vec(d, 1, 0) * sin(t) +
        ref_matrix_eig(d, 1) * ref_matrix_vec(d, 1, 1) * cos(t);

    z = origin[2] + ref_matrix_eig(d, 0) * ref_matrix_vec(d, 2, 0) * sin(t) +
        ref_matrix_eig(d, 1) * ref_matrix_vec(d, 2, 1) * cos(t);

    fprintf(f, " %.15e %.15e %.15e\n", x, y, z);
  }
  for (i = 0; i < n; i++) {
    t = dt * (REF_DBL)i;

    x = origin[0] + ref_matrix_eig(d, 0) * ref_matrix_vec(d, 0, 0) * sin(t) +
        ref_matrix_eig(d, 2) * ref_matrix_vec(d, 0, 2) * cos(t);

    y = origin[1] + ref_matrix_eig(d, 0) * ref_matrix_vec(d, 1, 0) * sin(t) +
        ref_matrix_eig(d, 2) * ref_matrix_vec(d, 1, 2) * cos(t);

    z = origin[2] + ref_matrix_eig(d, 0) * ref_matrix_vec(d, 2, 0) * sin(t) +
        ref_matrix_eig(d, 2) * ref_matrix_vec(d, 2, 2) * cos(t);

    fprintf(f, " %.15e %.15e %.15e\n", x, y, z);
  }
  for (i = 0; i < n; i++) {
    t = dt * (REF_DBL)i;

    x = origin[0] + ref_matrix_eig(d, 1) * ref_matrix_vec(d, 0, 1) * sin(t) +
        ref_matrix_eig(d, 2) * ref_matrix_vec(d, 0, 2) * cos(t);

    y = origin[1] + ref_matrix_eig(d, 1) * ref_matrix_vec(d, 1, 1) * sin(t) +
        ref_matrix_eig(d, 2) * ref_matrix_vec(d, 1, 2) * cos(t);

    z = origin[2] + ref_matrix_eig(d, 1) * ref_matrix_vec(d, 2, 1) * sin(t) +
        ref_matrix_eig(d, 2) * ref_matrix_vec(d, 2, 2) * cos(t);

    fprintf(f, " %.15e %.15e %.15e\n", x, y, z);
  }
  fprintf(f, "      ' />\n");
  fprintf(f, "      </IndexedLineSet>\n");

  return REF_SUCCESS;
}

REF_STATUS ref_html_free(REF_HTML ref_html) {
  if (NULL != (void *)ref_html_file(ref_html)) {
    fprintf(ref_html_file(ref_html), "    </shape></scene></x3d>\n");
    fprintf(ref_html_file(ref_html), "  </body>\n");
    fprintf(ref_html_file(ref_html), "</html>\n");
    fclose(ref_html_file(ref_html));
    ref_html_file(ref_html) = NULL;
  }
  ref_free(ref_html);
  return REF_SUCCESS;
}
