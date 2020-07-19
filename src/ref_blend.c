
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

#include "ref_blend.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_dict.h"
#include "ref_egads.h"
#include "ref_malloc.h"

#define ref_blend_geom(ref_blend) ((ref_blend)->geom)
#define ref_blend_grid(ref_blend) ((ref_blend)->grid)
#define ref_blend_displacement(ref_blend, ixyz, geom) \
  ((ref_blend)->displacement[(ixyz) + 3 * (geom)])

REF_STATUS ref_blend_create(REF_BLEND *ref_blend_ptr, REF_GRID ref_grid) {
  REF_BLEND ref_blend;
  REF_INT n;

  ref_malloc(*ref_blend_ptr, 1, REF_BLEND_STRUCT);

  ref_blend = *ref_blend_ptr;

  ref_blend_geom(ref_blend) = ref_grid_geom(ref_grid);

  RSS(ref_grid_deep_copy(&ref_blend_grid(ref_blend), ref_grid), "deep copy");
  n = ref_geom_max(ref_grid_geom(ref_blend_grid(ref_blend)));
  ref_malloc_init(ref_blend->displacement, 3 * n, REF_DBL, 0.0);

  return REF_SUCCESS;
}

REF_STATUS ref_blend_free(REF_BLEND ref_blend) {
  if (NULL == (void *)ref_blend) return REF_NULL;

  ref_free(ref_blend->displacement);
  ref_grid_free(ref_blend_grid(ref_blend));
  /* geom is a pointer to the orig */
  ref_free(ref_blend);

  return REF_SUCCESS;
}

REF_STATUS ref_blend_initialize(REF_BLEND ref_blend) {
  REF_GEOM ref_geom = ref_blend_geom(ref_blend);
  REF_DBL edge_xyz[3], face_xyz[3];
  REF_INT i, edge_geom, face_geom, node, item;

  /* also displace edge and face to geom nodes */

  each_ref_geom_edge(ref_geom, edge_geom) {
    RSS(ref_egads_eval(ref_geom, edge_geom, edge_xyz, NULL), "eval edge");
    node = ref_geom_node(ref_geom, edge_geom);
    each_ref_geom_having_node(ref_geom, node, item, face_geom) {
      if (REF_GEOM_FACE == ref_geom_type(ref_geom, node)) {
        RSS(ref_egads_eval(ref_geom, face_geom, face_xyz, NULL), "eval face");
        for (i = 0; i < 3; i++) {
          ref_blend_displacement(ref_blend, i, face_geom) =
              edge_xyz[i] - face_xyz[i];
        }
      }
    }
  }
  return REF_SUCCESS;
}

static REF_STATUS ref_blend_face_tec_zone(REF_BLEND ref_blend, REF_INT id,
                                          FILE *file) {
  REF_GRID ref_grid = ref_blend_grid(ref_blend);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_GEOM ref_geom = ref_blend_geom(ref_blend);
  REF_DICT ref_dict, ref_dict_jump, ref_dict_degen;
  REF_INT geom, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local, node;
  REF_INT nnode, nnode_sens0, nnode_degen, ntri;
  REF_INT sens;
  REF_DBL *uv, param[2];
  REF_DBL xyz[3];

  RSS(ref_dict_create(&ref_dict), "create dict");
  RSS(ref_dict_create(&ref_dict_jump), "create dict");
  RSS(ref_dict_create(&ref_dict_degen), "create dict");

  each_ref_geom_face(ref_geom, geom) {
    node = ref_geom_node(ref_geom, geom);
    if (id == ref_geom_id(ref_geom, geom)) {
      if (0 == ref_geom_degen(ref_geom, geom)) {
        RSS(ref_dict_store(ref_dict, node, geom), "mark nodes");
        if (0 != ref_geom_jump(ref_geom, geom)) {
          RSS(ref_dict_store(ref_dict_jump, node, geom), "mark jump");
        }
      } else {
        each_ref_cell_having_node(ref_cell, node, item, cell) {
          RSS(ref_cell_nodes(ref_cell, cell, nodes), "nodes");
          if (id == nodes[3]) {
            RSS(ref_dict_store(ref_dict_degen, cell, node), "mark degen");
          }
        }
      }
    }
  }

  nnode_sens0 = ref_dict_n(ref_dict);
  nnode_degen = ref_dict_n(ref_dict) + ref_dict_n(ref_dict_jump);
  nnode = ref_dict_n(ref_dict) + ref_dict_n(ref_dict_jump) +
          ref_dict_n(ref_dict_degen);

  ntri = 0;
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[3]) {
      ntri++;
    }
  }

  /* skip degenerate */
  if (0 == nnode || 0 == ntri) {
    RSS(ref_dict_free(ref_dict_degen), "free degen");
    RSS(ref_dict_free(ref_dict_jump), "free jump");
    RSS(ref_dict_free(ref_dict), "free dict");
    return REF_SUCCESS;
  }

  fprintf(
      file,
      "zone t=\"face%d\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
      id, nnode, ntri, "point", "fetriangle");

  ref_malloc_init(uv, 2 * nnode, REF_DBL, -1.0);
  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[3]) {
      each_ref_cell_cell_node(ref_cell, node) {
        RSS(ref_geom_find(ref_geom, nodes[node], REF_GEOM_FACE, id, &geom),
            "find");
        RSS(ref_geom_cell_tuv(ref_geom, nodes[node], nodes, REF_GEOM_FACE,
                              param, &sens),
            "cell tuv");
        if (0 == ref_geom_degen(ref_geom, geom)) {
          if (0 == sens || 1 == sens) {
            RSS(ref_dict_location(ref_dict, nodes[node], &local), "localize");
          } else {
            RSS(ref_dict_location(ref_dict_jump, nodes[node], &local),
                "localize");
            local += nnode_sens0;
          }
        } else {
          RSS(ref_dict_location(ref_dict_degen, cell, &local), "localize");
          local += nnode_degen;
        }
        uv[0 + 2 * local] = param[0];
        uv[1 + 2 * local] = param[1];
      }
    }
  }

  each_ref_dict_key_value(ref_dict, item, node, geom) {
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);

    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], uv[0 + 2 * item], uv[1 + 2 * item],
            ref_blend_displacement(ref_blend, 0, geom),
            ref_blend_displacement(ref_blend, 1, geom),
            ref_blend_displacement(ref_blend, 2, geom));
  }
  each_ref_dict_key_value(ref_dict_jump, item, node, geom) {
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], uv[0 + 2 * (nnode_sens0 + item)],
            uv[1 + 2 * (nnode_sens0 + item)],
            ref_blend_displacement(ref_blend, 0, geom),
            ref_blend_displacement(ref_blend, 1, geom),
            ref_blend_displacement(ref_blend, 2, geom));
  }
  each_ref_dict_key_value(ref_dict_degen, item, cell, node) {
    xyz[0] = ref_node_xyz(ref_node, 0, node);
    xyz[1] = ref_node_xyz(ref_node, 1, node);
    xyz[2] = ref_node_xyz(ref_node, 2, node);
    fprintf(file, " %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", xyz[0],
            xyz[1], xyz[2], uv[0 + 2 * (nnode_degen + item)],
            uv[1 + 2 * (nnode_degen + item)],
            ref_blend_displacement(ref_blend, 0, geom),
            ref_blend_displacement(ref_blend, 1, geom),
            ref_blend_displacement(ref_blend, 2, geom));
  }
  ref_free(uv);

  each_ref_cell_valid_cell_with_nodes(ref_cell, cell, nodes) {
    if (id == nodes[3]) {
      each_ref_cell_cell_node(ref_cell, node) {
        RSS(ref_geom_find(ref_geom, nodes[node], REF_GEOM_FACE, id, &geom),
            "find");
        RSS(ref_geom_cell_tuv(ref_geom, nodes[node], nodes, REF_GEOM_FACE,
                              param, &sens),
            "cell tuv");
        if (0 == ref_geom_degen(ref_geom, geom)) {
          if (0 == sens || 1 == sens) {
            RSS(ref_dict_location(ref_dict, nodes[node], &local), "localize");
          } else {
            RSS(ref_dict_location(ref_dict_jump, nodes[node], &local),
                "localize");
            local += nnode_sens0;
          }
        } else {
          RSS(ref_dict_location(ref_dict_degen, cell, &local), "localize");
          local += nnode_degen;
        }
        fprintf(file, " %d", local + 1);
      }
      fprintf(file, "\n");
    }
  }

  RSS(ref_dict_free(ref_dict_degen), "free degen");
  RSS(ref_dict_free(ref_dict_jump), "free jump");
  RSS(ref_dict_free(ref_dict), "free dict");

  return REF_SUCCESS;
}

REF_STATUS ref_blend_tec(REF_BLEND ref_blend, const char *filename) {
  REF_GEOM ref_geom = ref_blend_geom(ref_blend);
  FILE *file;
  REF_INT geom, id, min_id, max_id;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"refine watertight cad displacement\"\n");
  fprintf(file,
          "variables = \"x\" \"y\" \"z\" \"u\" \"v\" \"dx\" \"dy\" \"dz\"\n");

  min_id = REF_INT_MAX;
  max_id = REF_INT_MIN;
  each_ref_geom_face(ref_geom, geom) {
    min_id = MIN(min_id, ref_geom_id(ref_geom, geom));
    max_id = MAX(max_id, ref_geom_id(ref_geom, geom));
  }

  for (id = min_id; id <= max_id; id++)
    RSS(ref_blend_face_tec_zone(ref_blend, id, file), "tec face");

  fclose(file);
  return REF_SUCCESS;
}
