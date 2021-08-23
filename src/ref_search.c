
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

#include "ref_search.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ref_malloc.h"
#include "ref_math.h"

#define MAX_NODE_LIST (100)

REF_STATUS ref_search_create(REF_SEARCH *ref_search_ptr, REF_INT n) {
  REF_SEARCH ref_search;

  ref_malloc(*ref_search_ptr, 1, REF_SEARCH_STRUCT);
  ref_search = (*ref_search_ptr);

  ref_search->d = 3;
  ref_search->n = n;
  ref_search->empty = 0;

  ref_malloc_init(ref_search->item, ref_search->n, REF_INT, REF_EMPTY);
  ref_malloc_init(ref_search->left, ref_search->n, REF_INT, REF_EMPTY);
  ref_malloc_init(ref_search->right, ref_search->n, REF_INT, REF_EMPTY);

  ref_malloc(ref_search->pos, ref_search->d * ref_search->n, REF_DBL);
  ref_malloc(ref_search->radius, ref_search->n, REF_DBL);
  ref_malloc_init(ref_search->children_ball, ref_search->n, REF_DBL, 0.0);

  return REF_SUCCESS;
}

REF_STATUS ref_search_free(REF_SEARCH ref_search) {
  if (NULL == (void *)ref_search) return REF_NULL;
  ref_free(ref_search->children_ball);
  ref_free(ref_search->radius);
  ref_free(ref_search->pos);
  ref_free(ref_search->right);
  ref_free(ref_search->left);
  ref_free(ref_search->item);
  ref_free(ref_search);
  return REF_SUCCESS;
}

static REF_STATUS ref_search_distance(REF_SEARCH ref_search, REF_INT a,
                                      REF_INT b, REF_DBL *distance) {
  REF_INT i;
  *distance = 0.0;
  for (i = 0; i < ref_search->d; i++)
    (*distance) += pow(ref_search->pos[i + ref_search->d * b] -
                           ref_search->pos[i + ref_search->d * a],
                       2);
  (*distance) = sqrt(*distance);
  return REF_SUCCESS;
}

static REF_STATUS ref_search_home(REF_SEARCH ref_search, REF_INT child,
                                  REF_INT parent) {
  REF_DBL child_distance;
  REF_DBL left_distance, right_distance;

  RUS(REF_EMPTY, child, "empty child");
  RUS(REF_EMPTY, parent, "empty parent");

  /* done, don't add self to children */
  if (child == parent) return REF_SUCCESS;

  RSS(ref_search_distance(ref_search, child, parent, &child_distance), "d");
  ref_search->children_ball[parent] =
      MAX(ref_search->children_ball[parent],
          child_distance + ref_search->radius[child]);

  if (REF_EMPTY == ref_search->left[parent]) {
    ref_search->left[parent] = child;
    return REF_SUCCESS;
  }

  if (REF_EMPTY == ref_search->right[parent]) {
    ref_search->right[parent] = child;
    return REF_SUCCESS;
  }

  RSS(ref_search_distance(ref_search, child, ref_search->left[parent],
                          &left_distance),
      "left dist");
  RSS(ref_search_distance(ref_search, child, ref_search->right[parent],
                          &right_distance),
      "right dist");

  if (left_distance < right_distance) {
    RSS(ref_search_home(ref_search, child, ref_search->left[parent]),
        "recursively add to left child");
  } else {
    RSS(ref_search_home(ref_search, child, ref_search->right[parent]),
        "recursively add to right child");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_search_insert(REF_SEARCH ref_search, REF_INT item,
                             REF_DBL *position, REF_DBL radius) {
  REF_INT i, location;
  if (ref_search->empty >= ref_search->n)
    RSS(REF_INCREASE_LIMIT, "need larger tree for more items");

  if (item < 0) RSS(REF_INVALID, "item can not be negative");

  location = ref_search->empty;
  (ref_search->empty)++;

  ref_search->item[location] = item;
  for (i = 0; i < ref_search->d; i++)
    ref_search->pos[i + ref_search->d * location] = position[i];
  ref_search->radius[location] = radius;

  RSS(ref_search_home(ref_search, location, 0), "top level home");

  return REF_SUCCESS;
}

static REF_STATUS ref_search_gather(REF_SEARCH ref_search, REF_LIST ref_list,
                                    REF_INT parent, REF_DBL *position,
                                    REF_DBL radius) {
  REF_INT i;
  REF_DBL distance;

  if (0 == ref_search->n) return REF_SUCCESS;  /* tree empty */
  if (REF_EMPTY == parent) return REF_SUCCESS; /* finished traversing */
  RAB(0 <= parent && parent < ref_search->n, "parent invalid",
      { printf("%d n %d parent\n", ref_search->n, parent); })
  /* finished traversing */
  if (REF_EMPTY == ref_search->item[parent]) return REF_SUCCESS;

  distance = 0.0;
  for (i = 0; i < ref_search->d; i++)
    distance +=
        pow(position[i] - ref_search->pos[i + ref_search->d * parent], 2);
  distance = sqrt(distance);

  /* if the distance between me and the target are less than combined radii */
  if (distance <= ref_search->radius[parent] + radius) {
    RSS(ref_list_push(ref_list, ref_search->item[parent]), "add item");
  }

  /* if the distance between me and the target are less than children
   * children_ball includes child radii, so only subtract target radius */
  if (distance - radius <= ref_search->children_ball[parent]) {
    RSS(ref_search_gather(ref_search, ref_list, ref_search->left[parent],
                          position, radius),
        "gthr");
    RSS(ref_search_gather(ref_search, ref_list, ref_search->right[parent],
                          position, radius),
        "gthr");
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_search_gather_tri(REF_SEARCH ref_search, REF_DBL *xyz,
                                        REF_INT parent, REF_DBL *position,
                                        REF_DBL *distance) {
  REF_INT i;
  REF_DBL dist;

  if (0 == ref_search->n) return REF_SUCCESS;  /* tree empty */
  if (REF_EMPTY == parent) return REF_SUCCESS; /* finished traversing */
  RAB(0 <= parent && parent < ref_search->n, "parent invalid",
      { printf("%d n %d parent\n", ref_search->n, parent); })
  /* finished traversing */
  if (REF_EMPTY == ref_search->item[parent]) return REF_SUCCESS;

  dist = 0.0;
  for (i = 0; i < ref_search->d; i++)
    dist += pow(position[i] - ref_search->pos[i + ref_search->d * parent], 2);
  dist = sqrt(dist);

  /* if the distance between me and the target are less than combined radii */
  if (dist - ref_search->radius[parent] <= *distance) {
    REF_INT tri = ref_search->item[parent];
    REF_DBL tri_dist;
    RSS(ref_search_distance3(&(xyz[0 + 9 * tri]), &(xyz[3 + 9 * tri]),
                             &(xyz[6 + 9 * tri]), position, &tri_dist),
        "tri dist");
    *distance = MIN(*distance, tri_dist);
  }

  /* if the distance between me and the target are less than children
   * children_ball includes child radii, so only subtract target radius */
  if (*distance >= dist - ref_search->children_ball[parent]) {
    RSS(ref_search_gather_tri(ref_search, xyz, ref_search->left[parent],
                              position, distance),
        "gthr");
    RSS(ref_search_gather_tri(ref_search, xyz, ref_search->right[parent],
                              position, distance),
        "gthr");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_search_touching(REF_SEARCH ref_search, REF_LIST ref_list,
                               REF_DBL *position, REF_DBL radius) {
  RSS(ref_search_gather(ref_search, ref_list, 0, position, radius), "gthr");
  return REF_SUCCESS;
}

static REF_STATUS ref_search_trim(REF_SEARCH ref_search, REF_INT parent,
                                  REF_DBL *position, REF_DBL *trim_radius) {
  REF_INT i;
  REF_DBL distance;

  if (REF_EMPTY == parent) return REF_SUCCESS;
  if (parent >= ref_search->n) return REF_SUCCESS;
  if (REF_EMPTY == ref_search->item[parent]) return REF_SUCCESS;

  distance = 0.0;
  for (i = 0; i < ref_search->d; i++)
    distance +=
        pow(position[i] - ref_search->pos[i + ref_search->d * parent], 2);
  distance = sqrt(distance);

  if (distance + ref_search->radius[parent] < *trim_radius) {
    *trim_radius = distance + ref_search->radius[parent];
  }

  /* if the trim_distance is larger than the distance between me and the target
   * minus the children_ball look for better */
  if (*trim_radius > distance - ref_search->children_ball[parent]) {
    RSS(ref_search_trim(ref_search, ref_search->left[parent], position,
                        trim_radius),
        "gthr");
    RSS(ref_search_trim(ref_search, ref_search->right[parent], position,
                        trim_radius),
        "gthr");
  }

  return REF_SUCCESS;
}

REF_STATUS ref_search_trim_radius(REF_SEARCH ref_search, REF_DBL *position,
                                  REF_DBL *trim_radius) {
  REF_INT parent;
  parent = 0;
  *trim_radius = REF_DBL_MAX;
  RSS(ref_search_trim(ref_search, parent, position, trim_radius), "trim");
  return REF_SUCCESS;
}

REF_STATUS ref_search_nearest_candidates(REF_SEARCH ref_search,
                                         REF_LIST ref_list, REF_DBL *position) {
  REF_DBL trim_radius;
  RSS(ref_search_trim_radius(ref_search, position, &trim_radius), "scope");
  RSS(ref_search_touching(ref_search, ref_list, position, trim_radius),
      "touches");
  return REF_SUCCESS;
}

REF_STATUS ref_search_nearest_candidates_closer_than(REF_SEARCH ref_search,
                                                     REF_LIST ref_list,
                                                     REF_DBL *position,
                                                     REF_DBL distance) {
  REF_DBL trim_radius;
  REF_INT parent;
  parent = 0;
  trim_radius = distance;
  RSS(ref_search_trim(ref_search, parent, position, &trim_radius), "trim");
  RSS(ref_search_touching(ref_search, ref_list, position, trim_radius),
      "touches");
  return REF_SUCCESS;
}
REF_STATUS ref_search_nearest_tri(REF_SEARCH ref_search, REF_DBL *xyz,
                                  REF_DBL *position, REF_DBL *distance) {
  REF_INT parent;
  parent = 0;
  RSS(ref_search_gather_tri(ref_search, xyz, parent, position, distance),
      "touches");
  return REF_SUCCESS;
}

REF_STATUS ref_search_selection(REF_MPI ref_mpi, REF_INT n, REF_DBL *elements,
                                REF_LONG position, REF_DBL *value) {
  REF_INT i, bisection;
  REF_LONG low_pos, high_pos, count;
  REF_DBL low_val, high_val, temp, mid_val;
  low_pos = 0;
  high_pos = (REF_LONG)n;
  RSS(ref_mpi_allsum(ref_mpi, &high_pos, 1, REF_LONG_TYPE), "high_pos");
  high_pos--;
  low_val = REF_DBL_MAX;
  high_val = REF_DBL_MIN;
  for (i = 0; i < n; i++) {
    low_val = MIN(low_val, elements[i]);
    high_val = MAX(high_val, elements[i]);
  }
  temp = low_val;
  RSS(ref_mpi_min(ref_mpi, &temp, &low_val, REF_DBL_TYPE), "min");
  RSS(ref_mpi_bcast(ref_mpi, &low_val, 1, REF_DBL_TYPE), "bcast");
  temp = high_val;
  RSS(ref_mpi_max(ref_mpi, &temp, &high_val, REF_DBL_TYPE), "max");
  RSS(ref_mpi_bcast(ref_mpi, &high_val, 1, REF_DBL_TYPE), "bcast");

  if (position <= low_pos) {
    *value = low_val;
    return REF_SUCCESS;
  }

  if (position >= high_pos) {
    *value = high_val;
    return REF_SUCCESS;
  }

  mid_val = 0.5 * (low_val + high_val); /* ensure initialized */
  for (bisection = 0; bisection < 40; bisection++) {
    mid_val = 0.5 * (low_val + high_val);
    count = 0;
    for (i = 0; i < n; i++) {
      if (elements[i] <= mid_val) count++;
    }

    RSS(ref_mpi_allsum(ref_mpi, &count, 1, REF_LONG_TYPE), "bcast");
    /* printf("pos  %ld %ld %ld val %f %f %f\n",
               low_pos, count, high_pos, low_val,
               mid_val, high_val);*/
    if (count - 1 < position) {
      low_val = mid_val;
    } else {
      high_val = mid_val;
    }
  }
  *value = mid_val;
  return REF_SUCCESS;
}

REF_STATUS ref_search_distance2(REF_DBL *xyz0, REF_DBL *xyz1, REF_DBL *xyz,
                                REF_DBL *distance) {
  REF_DBL dl[3], dxyz[3], len2, proj2, t;
  dl[0] = xyz1[0] - xyz0[0];
  dl[1] = xyz1[1] - xyz0[1];
  dl[2] = xyz1[2] - xyz0[2];
  dxyz[0] = xyz[0] - xyz0[0];
  dxyz[1] = xyz[1] - xyz0[1];
  dxyz[2] = xyz[2] - xyz0[2];
  len2 = ref_math_dot(dl, dl);
  proj2 = ref_math_dot(dxyz, dl);
  if (ref_math_divisible(proj2, len2)) {
    t = proj2 / len2;
    t = MAX(t, 0.0);
    t = MIN(t, 1.0);
    dxyz[0] = xyz[0] - (xyz0[0] + t * dl[0]);
    dxyz[1] = xyz[1] - (xyz0[1] + t * dl[1]);
    dxyz[2] = xyz[2] - (xyz0[2] + t * dl[2]);
    *distance = sqrt(ref_math_dot(dxyz, dxyz));
  } else { /* length zero, either endpoint */
    *distance = sqrt(ref_math_dot(dxyz, dxyz));
  }

  return REF_SUCCESS;
}

static REF_STATUS ref_search_xyz_normal(REF_DBL *xyz0, REF_DBL *xyz1,
                                        REF_DBL *xyz2, REF_DBL *normal) {
  REF_DBL edge10[3], edge20[3];

  edge10[0] = xyz1[0] - xyz0[0];
  edge10[1] = xyz1[1] - xyz0[1];
  edge10[2] = xyz1[2] - xyz0[2];

  edge20[0] = xyz2[0] - xyz0[0];
  edge20[1] = xyz2[1] - xyz0[1];
  edge20[2] = xyz2[2] - xyz0[2];

  ref_math_cross_product(edge10, edge20, normal);

  return REF_SUCCESS;
}

REF_STATUS ref_search_distance3(REF_DBL *xyz0, REF_DBL *xyz1, REF_DBL *xyz2,
                                REF_DBL *xyz, REF_DBL *distance) {
  REF_DBL dist;
  REF_DBL bary[3], total, total_normal[3], normal[3];
  REF_DBL dxyz[3];
  REF_DBL xyzp[3];

  RSS(ref_search_xyz_normal(xyz0, xyz1, xyz2, total_normal), "n0");

  /* projects query point to triangle plane */
  xyzp[0] = xyz[0] - xyz0[0];
  xyzp[1] = xyz[1] - xyz0[1];
  xyzp[2] = xyz[2] - xyz0[2];
  total = ref_math_dot(xyzp, total_normal);
  xyzp[0] -= total_normal[0] * total;
  xyzp[1] -= total_normal[1] * total;
  xyzp[2] -= total_normal[2] * total;
  xyzp[0] += xyz0[0];
  xyzp[1] += xyz0[1];
  xyzp[2] += xyz0[2];

  RSS(ref_search_xyz_normal(xyzp, xyz1, xyz2, normal), "n0");
  bary[0] = ref_math_dot(normal, total_normal);
  RSS(ref_search_xyz_normal(xyz0, xyzp, xyz2, normal), "n1");
  bary[1] = ref_math_dot(normal, total_normal);
  RSS(ref_search_xyz_normal(xyz0, xyz1, xyzp, normal), "n2");
  bary[2] = ref_math_dot(normal, total_normal);

  total = bary[0] + bary[1] + bary[2];

  if (ref_math_divisible(bary[0], total) &&
      ref_math_divisible(bary[1], total) &&
      ref_math_divisible(bary[2], total)) {
    bary[0] /= total;
    bary[1] /= total;
    bary[2] /= total;
    if (bary[0] >= 0.0 && bary[1] >= 0.0 && bary[2] >= 0.0) {
      dxyz[0] =
          bary[0] * xyz0[0] + bary[1] * xyz1[0] + bary[2] * xyz2[0] - xyz[0];
      dxyz[1] =
          bary[0] * xyz0[1] + bary[1] * xyz1[1] + bary[2] * xyz2[1] - xyz[1];
      dxyz[2] =
          bary[0] * xyz0[2] + bary[1] * xyz1[2] + bary[2] * xyz2[2] - xyz[2];
      *distance = sqrt(ref_math_dot(dxyz, dxyz));
      return REF_SUCCESS;
    }
  }

  RSS(ref_search_distance2(xyz0, xyz1, xyz, &dist), "e01");
  *distance = dist;
  RSS(ref_search_distance2(xyz1, xyz2, xyz, &dist), "e12");
  *distance = MIN(*distance, dist);
  RSS(ref_search_distance2(xyz2, xyz0, xyz, &dist), "e20");
  *distance = MIN(*distance, dist);

  return REF_SUCCESS;
}
/* Ericson Real Time Collision Detection p141 */
REF_STATUS ref_search_dist3(REF_DBL *a, REF_DBL *b, REF_DBL *c, REF_DBL *p,
                            REF_DBL *distance) {
  REF_DBL ab[3], ac[3], ap[3];
  REF_DBL d1, d2;
  REF_DBL bp[3];
  REF_DBL d3, d4;
  REF_DBL vc, v;
  REF_DBL proj[3];
  REF_DBL cp[3];
  REF_DBL d5, d6;
  REF_DBL vb, va;
  REF_DBL denom, w;

  /* Check if P in vertex region outside A */
  ab[0] = b[0] - a[0];
  ab[1] = b[1] - a[1];
  ab[2] = b[2] - a[2];

  ac[0] = c[0] - a[0];
  ac[1] = c[1] - a[1];
  ac[2] = c[2] - a[2];

  ap[0] = p[0] - a[0];
  ap[1] = p[1] - a[1];
  ap[2] = p[2] - a[2];

  d1 = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2];
  d2 = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2];
  if (d1 <= 0.0 && d2 <= 0.0) {
    *distance =
        sqrt((p[0] - a[0]) * (p[0] - a[0]) + (p[1] - a[1]) * (p[1] - a[1]) +
             (p[2] - a[2]) * (p[2] - a[2]));
    return REF_SUCCESS;
  }

  /* Check if P in vertex region outside B */
  bp[0] = p[0] - b[0];
  bp[1] = p[1] - b[1];
  bp[2] = p[2] - b[2];

  d3 = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2];
  d4 = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2];
  if (d3 >= 0.0 && d4 <= d3) {
    *distance =
        sqrt((p[0] - b[0]) * (p[0] - b[0]) + (p[1] - b[1]) * (p[1] - b[1]) +
             (p[2] - b[2]) * (p[2] - b[2]));
    return REF_SUCCESS;
  }

  /* Check if P in edge region of AB, if so return projection of P onto AB */
  vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    RAS(ref_math_divisible(d1, (d1 - d3)), "div zero d1/(d1-d3)");
    v = d1 / (d1 - d3);
    proj[0] = a[0] + v * ab[0];
    proj[1] = a[1] + v * ab[1];
    proj[2] = a[2] + v * ab[2];
    *distance = sqrt((p[0] - proj[0]) * (p[0] - proj[0]) +
                     (p[1] - proj[1]) * (p[1] - proj[1]) +
                     (p[2] - proj[2]) * (p[2] - proj[2]));
    return REF_SUCCESS;
  }

  /* Check if P in vertex region outside C */
  cp[0] = p[0] - c[0];
  cp[1] = p[1] - c[1];
  cp[2] = p[2] - c[2];
  d5 = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2];
  d6 = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2];
  if (d6 >= 0.0 && d5 <= d6) {
    *distance =
        sqrt((p[0] - c[0]) * (p[0] - c[0]) + (p[1] - c[1]) * (p[1] - c[1]) +
             (p[2] - c[2]) * (p[2] - c[2]));
    return REF_SUCCESS;
  }

  /* Check if P in edge region of AC, if so return projection of P onto AC */
  vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    RAS(ref_math_divisible(d1, (d1 - d3)), "div zero d1/(d1-d3)");
    v = d2 / (d2 - d6);
    proj[0] = a[0] + v * ac[0];
    proj[1] = a[1] + v * ac[1];
    proj[2] = a[2] + v * ac[2];
    *distance = sqrt((p[0] - proj[0]) * (p[0] - proj[0]) +
                     (p[1] - proj[1]) * (p[1] - proj[1]) +
                     (p[2] - proj[2]) * (p[2] - proj[2]));
    return REF_SUCCESS;
  }

  /* Check if P in edge region of BC, if so return projection of P onto BC */
  va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    RAS(ref_math_divisible((d4 - d3), ((d4 - d3) + (d5 - d6))),
        "div zero (d4 - d3) / ((d4 - d3) + (d5 - d6))");
    v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    proj[0] = b[0] + v * (c[0] - b[0]);
    proj[1] = b[1] + v * (c[1] - b[1]);
    proj[2] = b[2] + v * (c[2] - b[2]);
    *distance = sqrt((p[0] - proj[0]) * (p[0] - proj[0]) +
                     (p[1] - proj[1]) * (p[1] - proj[1]) +
                     (p[2] - proj[2]) * (p[2] - proj[2]));
    return REF_SUCCESS;
  }

  RAS(ref_math_divisible(1.0, (va + vb + vc)), "div zero 1.0 / (va + vb + vc)");
  denom = 1.0 / (va + vb + vc);
  v = vb * denom;
  w = vc * denom;
  proj[0] = a[0] + v * ab[0] + w * ac[0];
  proj[1] = a[1] + v * ab[1] + w * ac[1];
  proj[2] = a[2] + v * ab[2] + w * ac[2];
  *distance = sqrt((p[0] - proj[0]) * (p[0] - proj[0]) +
                   (p[1] - proj[1]) * (p[1] - proj[1]) +
                   (p[2] - proj[2]) * (p[2] - proj[2]));

  return REF_SUCCESS;
}
