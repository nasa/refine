
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

REF_STATUS ref_search_distance(REF_SEARCH ref_search, REF_INT a, REF_INT b,
                               REF_DBL *distance) {
  REF_INT i;
  *distance = 0.0;
  for (i = 0; i < ref_search->d; i++)
    (*distance) += pow(ref_search->pos[i + ref_search->d * b] -
                           ref_search->pos[i + ref_search->d * a],
                       2);
  (*distance) = sqrt(*distance);
  return REF_SUCCESS;
}

REF_STATUS ref_search_home(REF_SEARCH ref_search, REF_INT child,
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

REF_STATUS ref_search_gather(REF_SEARCH ref_search, REF_LIST ref_list,
                             REF_INT parent, REF_DBL *position,
                             REF_DBL radius) {
  REF_INT i;
  REF_DBL distance;

  if (REF_EMPTY == parent) return REF_SUCCESS;
  RAS(0 <= parent && parent < ref_search->n, "parent invalid")
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
