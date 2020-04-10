
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

#include "ref_adj.h"
#include "ref_args.h"
#include "ref_cell.h"
#include "ref_dict.h"
#include "ref_edge.h"
#include "ref_export.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_inflate.h"
#include "ref_list.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_node.h"
#include "ref_part.h"
#include "ref_sort.h"

int main(int argc, char *argv[]) {
  REF_MPI ref_mpi;
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_DICT faceids;
  REF_INT arg, faceid, nlayers;
  REF_DBL first_thickness, total_thickness, mach;
  REF_DBL rate, total;
  REF_INT layer;
  REF_DBL thickness, xshift, mach_angle_rad;
  REF_INT aoa_pos;
  REF_DBL alpha_deg = 0;
  REF_DBL alpha_rad = 0;
  REF_INT origin_pos;
  REF_INT rotate_pos;
  REF_DBL rotate_deg = 0;
  REF_DBL rotate_rad = 0;
  REF_INT scale_pos;
  REF_DBL scale = 0;
  REF_INT pos;
  REF_INT node;
  REF_DBL x, z;
  REF_BOOL extrude_radially = REF_FALSE;
  REF_DBL origin[3];
  REF_INT last_face_arg;
  REF_INT mapbc_pos;
  REF_INT bc_type;
  char *mapbc_file_name, *family_name;

  if (7 > argc) {
    printf(
        "usage: \n %s input.grid nlayers first_thickness total_thickness mach "
        "faceid  [faceid...]\n",
        argv[0]);
    printf("       [--mapbc usm3d_format.mapbc family_name bc_type]\n");
    printf("       [--aoa angle_of_attack_in_degrees]\n");
    printf("       [--rotate angle_in_degrees] (applied before inflation)\n");
    printf("       [--shift dx dy dz] (applied before inflation)\n");
    printf("       [--origin ox oy oz]\n");
    printf("       [--scale factor] (applied after inflation)\n");
    printf("  when first_thickness <= 0, it is set to a uniform grid,\n");
    printf("    first_thickness = total_thickness/nlayers\n");
    printf("  when nlayers < 0, extrude radially\n");
    printf("    (--aoa option only available for radial extrusion)\n");
    return 1;
  }

  RSS(ref_mpi_start(argc, argv), "start");
  RSS(ref_mpi_create(&ref_mpi), "create");
  ref_mpi_stopwatch_start(ref_mpi);

  if (ref_mpi_para(ref_mpi)) {
    if (ref_mpi_once(ref_mpi)) printf(" part %s\n", argv[1]);
    RSS(ref_part_by_extension(&ref_grid, ref_mpi, argv[1]), "part");
    ref_mpi_stopwatch_stop(ref_mpi, "part grid");
    RSS(ref_migrate_to_balance(ref_grid), "balance");
    ref_mpi_stopwatch_stop(ref_mpi, "balance grid");
    RSS(ref_grid_pack(ref_grid), "pack");
    ref_mpi_stopwatch_stop(ref_mpi, "pack grid");
  } else {
    if (ref_mpi_once(ref_mpi)) printf(" import %s\n", argv[1]);
    RSS(ref_import_by_extension(&ref_grid, ref_mpi, argv[1]), "import");
  }
  ref_mpi_stopwatch_stop(ref_mpi, "read grid");

  nlayers = atoi(argv[2]);
  if (nlayers < 0) {
    nlayers = ABS(nlayers);
    extrude_radially = REF_TRUE;
  }
  first_thickness = atof(argv[3]);
  total_thickness = atof(argv[4]);
  mach = atof(argv[5]);

  last_face_arg = argc;

  aoa_pos = REF_EMPTY;
  RXS(ref_args_find(argc, argv, "--aoa", &aoa_pos), REF_NOT_FOUND,
      "aoa search");

  if (REF_EMPTY != aoa_pos) {
    if (aoa_pos >= argc) THROW("--aoa requires a value");
    if (!extrude_radially)
      THROW("--aoa requires radial extrusion, nlayers < 0");
    alpha_deg = atof(argv[aoa_pos + 1]);
    alpha_rad = ref_math_in_radians(alpha_deg);
    if (ref_mpi_once(ref_mpi)) printf(" --aoa %f deg\n", alpha_deg);
    last_face_arg = MIN(last_face_arg, aoa_pos);
  }

  origin_pos = REF_EMPTY;
  RXS(ref_args_find(argc, argv, "--origin", &origin_pos), REF_NOT_FOUND,
      "origin search");

  if (REF_EMPTY != origin_pos) {
    if (origin_pos >= argc - 3) THROW("--origin requires a value");
    origin[0] = atof(argv[origin_pos + 1]);
    origin[1] = atof(argv[origin_pos + 2]);
    origin[2] = atof(argv[origin_pos + 3]);
    if (ref_mpi_once(ref_mpi))
      printf(" --origin %f %f %f\n", origin[0], origin[1], origin[2]);
    last_face_arg = MIN(last_face_arg, origin_pos);
  }

  pos = REF_EMPTY;
  RXS(ref_args_find(argc, argv, "--shift", &pos), REF_NOT_FOUND,
      "shift search");

  if (REF_EMPTY != pos) {
    REF_DBL dx, dy, dz;
    if (pos >= argc - 3) THROW("--shift requires three values");
    dx = atof(argv[pos + 1]);
    dy = atof(argv[pos + 2]);
    dz = atof(argv[pos + 3]);
    if (ref_mpi_once(ref_mpi)) printf(" --shift %f %f %f\n", dx, dy, dz);
    last_face_arg = MIN(last_face_arg, pos);

    ref_node = ref_grid_node(ref_grid);
    each_ref_node_valid_node(ref_node, node) {
      ref_node_xyz(ref_node, 0, node) += dx;
      ref_node_xyz(ref_node, 1, node) += dy;
      ref_node_xyz(ref_node, 2, node) += dz;
    }
  }

  rotate_pos = REF_EMPTY;
  RXS(ref_args_find(argc, argv, "--rotate", &rotate_pos), REF_NOT_FOUND,
      "rotate search");

  if (REF_EMPTY != rotate_pos) {
    if (rotate_pos >= argc - 1) THROW("--rotate requires a value");
    rotate_deg = atof(argv[rotate_pos + 1]);
    rotate_rad = ref_math_in_radians(rotate_deg);
    if (ref_mpi_once(ref_mpi))
      printf(" --rotate %f deg (%f rad)\n", rotate_deg, rotate_rad);
    last_face_arg = MIN(last_face_arg, rotate_pos);

    ref_node = ref_grid_node(ref_grid);
    each_ref_node_valid_node(ref_node, node) {
      x = ref_node_xyz(ref_node, 0, node);
      z = ref_node_xyz(ref_node, 2, node);
      ref_node_xyz(ref_node, 0, node) =
          x * cos(rotate_rad) - z * sin(rotate_rad);
      ref_node_xyz(ref_node, 2, node) =
          x * sin(rotate_rad) + z * cos(rotate_rad);
    }
  }

  scale_pos = REF_EMPTY;
  RXS(ref_args_find(argc, argv, "--scale", &scale_pos), REF_NOT_FOUND,
      "scale search");

  if (REF_EMPTY != scale_pos) {
    if (scale_pos >= argc - 1) THROW("--scale requires a value");
    scale = atof(argv[scale_pos + 1]);
    if (ref_mpi_once(ref_mpi)) printf(" --scale %f\n", scale);
    last_face_arg = MIN(last_face_arg, scale_pos);
  }

  mapbc_pos = REF_EMPTY;
  RXS(ref_args_find(argc, argv, "--mapbc", &mapbc_pos), REF_NOT_FOUND,
      "mapbc search");

  RSS(ref_dict_create(&faceids), "create");

  if (REF_EMPTY != mapbc_pos) {
    if (mapbc_pos >= argc - 3) THROW("--mapbc requires three values");
    mapbc_file_name = argv[mapbc_pos + 1];
    family_name = argv[mapbc_pos + 2];
    bc_type = atoi(argv[mapbc_pos + 3]);
    if (ref_mpi_once(ref_mpi))
      printf(" --mapbc %s %s %d\n", mapbc_file_name, family_name, bc_type);
    last_face_arg = MIN(last_face_arg, mapbc_pos);
    RSS(ref_inflate_read_usm3d_mapbc(faceids, mapbc_file_name, family_name,
                                     bc_type),
        "faceids from mapbc");
  }

  if (ref_mpi_once(ref_mpi)) printf("faceids\n");
  for (arg = 6; arg < last_face_arg; arg++) {
    faceid = atoi(argv[arg]);
    RSS(ref_dict_store(faceids, faceid, REF_EMPTY), "store");
  }
  if (ref_mpi_once(ref_mpi))
    RSS(ref_dict_inspect_keys(faceids), "faceids dict inspect");

  if (first_thickness <= 0.0) {
    first_thickness = total_thickness / (REF_DBL)nlayers;
    rate = 1.0;
  } else {
    RSS(ref_inflate_rate(nlayers, first_thickness, total_thickness, &rate),
        "compute rate");
  }

  mach_angle_rad = asin(1 / mach);

  if (ref_mpi_once(ref_mpi)) {
    printf("inflating %d faces\n", ref_dict_n(faceids));
    printf("mach %f mach angle %f rad %f deg\n", mach, mach_angle_rad,
           ref_math_in_degrees(mach_angle_rad));
    printf("first thickness %f\n", first_thickness);
    printf("total thickness %f\n", total_thickness);
    printf("rate %f\n", rate);
    printf("layers %d\n", nlayers);
    printf("extrusion %d\n", extrude_radially);
  }

  if (REF_EMPTY == origin_pos)
    RSS(ref_inflate_origin(ref_grid, faceids, origin), "orig");

  total = 0.0;
  for (layer = 0; layer < nlayers; layer++) {
    thickness = first_thickness * pow(rate, layer);
    total = total + thickness;
    xshift = thickness / tan(mach_angle_rad);
    if (extrude_radially) {
      RSS(ref_inflate_radially(ref_grid, faceids, origin, thickness,
                               mach_angle_rad, alpha_rad),
          "inflate");
    } else {
      RSS(ref_inflate_face(ref_grid, faceids, origin, thickness, xshift),
          "inflate");
    }
    if (ref_mpi_once(ref_mpi))
      printf("layer%5d of%5d thickness %10.3e total %10.3e " REF_GLOB_FMT
             " nodes\n",
             layer + 1, nlayers, thickness, total,
             ref_node_n_global(ref_grid_node(ref_grid)));
  }

  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "inflate");

  if (ref_mpi_once(ref_mpi)) {
    printf("inflating %d faces\n", ref_dict_n(faceids));
    printf("mach %f mach angle %f rad %f deg\n", mach, mach_angle_rad,
           ref_math_in_degrees(mach_angle_rad));
    printf("first thickness %f\n", first_thickness);
    printf("total thickness %f\n", total_thickness);
    printf("rate %f\n", rate);
    printf("layers %d\n", nlayers);
    printf("extrusion %d\n", extrude_radially);
  }

  if (REF_EMPTY != scale_pos) {
    if (ref_mpi_once(ref_mpi))
      printf("scale grid after inflation by %f\n", scale);
    ref_node = ref_grid_node(ref_grid);
    each_ref_node_valid_node(ref_node, node) {
      ref_node_xyz(ref_node, 0, node) *= scale;
      ref_node_xyz(ref_node, 1, node) *= scale;
      ref_node_xyz(ref_node, 2, node) *= scale;
    }
  }

  RSS(ref_gather_by_extension(ref_grid, "inflated.b8.ugrid"), "b8");
  ref_mpi_stopwatch_stop(ref_grid_mpi(ref_grid), "export");

  RSS(ref_dict_free(faceids), "free");
  RSS(ref_grid_free(ref_grid), "free");

  ref_mpi_stopwatch_stop(ref_mpi, "done.");
  RSS(ref_mpi_free(ref_mpi), "free");
  RSS(ref_mpi_stop(), "stop");

  return 0;
}
