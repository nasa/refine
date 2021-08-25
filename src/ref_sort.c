
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

#include "ref_sort.h"

#include <stdio.h>
#include <stdlib.h>

#include "ref_malloc.h"

REF_STATUS ref_sort_insertion_int(REF_INT n, REF_INT *original,
                                  REF_INT *sorted) {
  REF_INT i, j, smallest, temp;

  for (i = 0; i < n; i++) sorted[i] = original[i];

  for (i = 0; i < n; i++) {
    smallest = i;
    for (j = i + 1; j < n; j++) {
      if (sorted[j] < sorted[smallest]) smallest = j;
    }
    temp = sorted[i];
    sorted[i] = sorted[smallest];
    sorted[smallest] = temp;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_sort_heap_int(REF_INT n, REF_INT *original,
                             REF_INT *sorted_index) {
  REF_INT i, j, l, ir, indxt, q;

  for (i = 0; i < n; i++) sorted_index[i] = i;

  if (n < 2) return REF_SUCCESS;

  l = (n >> 1) + 1;
  ir = n - 1;
  for (;;) {
    if (l > 1) {
      l--;
      indxt = sorted_index[l - 1];
      q = original[indxt];
    } else {
      indxt = sorted_index[ir];
      q = original[indxt];
      sorted_index[ir] = sorted_index[0];
      if (--ir == 0) {
        sorted_index[0] = indxt;
        break;
      }
    }
    i = l - 1;
    j = l + i;

    while (j <= ir) {
      if (j < ir) {
        if (original[sorted_index[j]] < original[sorted_index[j + 1]]) j++;
      }
      if (q < original[sorted_index[j]]) {
        sorted_index[i] = sorted_index[j];
        i = j;

        j++;
        j <<= 1;
        j--;

      } else
        break;
    }
    sorted_index[i] = indxt;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_sort_heap_glob(REF_INT n, REF_GLOB *original,
                              REF_INT *sorted_index) {
  REF_GLOB q;
  REF_INT i, j, l, ir, indxt;

  for (i = 0; i < n; i++) sorted_index[i] = i;

  if (n < 2) return REF_SUCCESS;

  l = (n >> 1) + 1;
  ir = n - 1;
  for (;;) {
    if (l > 1) {
      l--;
      indxt = sorted_index[l - 1];
      q = original[indxt];
    } else {
      indxt = sorted_index[ir];
      q = original[indxt];
      sorted_index[ir] = sorted_index[0];
      if (--ir == 0) {
        sorted_index[0] = indxt;
        break;
      }
    }
    i = l - 1;
    j = l + i;

    while (j <= ir) {
      if (j < ir) {
        if (original[sorted_index[j]] < original[sorted_index[j + 1]]) j++;
      }
      if (q < original[sorted_index[j]]) {
        sorted_index[i] = sorted_index[j];
        i = j;

        j++;
        j <<= 1;
        j--;

      } else
        break;
    }
    sorted_index[i] = indxt;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_sort_heap_dbl(REF_INT n, REF_DBL *original,
                             REF_INT *sorted_index) {
  REF_INT i, j, l, ir, indxt;
  REF_DBL q;

  for (i = 0; i < n; i++) sorted_index[i] = i;

  if (n < 2) return REF_SUCCESS;

  l = (n >> 1) + 1;
  ir = n - 1;
  for (;;) {
    if (l > 1) {
      l--;
      indxt = sorted_index[l - 1];
      q = original[indxt];
    } else {
      indxt = sorted_index[ir];
      q = original[indxt];
      sorted_index[ir] = sorted_index[0];
      if (--ir == 0) {
        sorted_index[0] = indxt;
        break;
      }
    }
    i = l - 1;
    j = l + i;

    while (j <= ir) {
      if (j < ir) {
        if (original[sorted_index[j]] < original[sorted_index[j + 1]]) j++;
      }
      if (q < original[sorted_index[j]]) {
        sorted_index[i] = sorted_index[j];
        i = j;

        j++;
        j <<= 1;
        j--;

      } else
        break;
    }
    sorted_index[i] = indxt;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_sort_in_place_glob(REF_INT n, REF_GLOB *sorts) {
  REF_INT i;
  REF_INT *order;
  REF_GLOB *sorted;

  /* see if it is too short to require sorting */
  if (2 > n) return REF_SUCCESS;
  ref_malloc(order, n, REF_INT);
  ref_malloc(sorted, n, REF_GLOB);
  RSS(ref_sort_heap_glob(n, sorts, order), "heap");
  for (i = 0; i < n; i++) {
    sorted[i] = sorts[order[i]];
  }
  for (i = 0; i < n; i++) {
    sorts[i] = sorted[i];
  }
  ref_free(sorted);
  ref_free(order);
  return REF_SUCCESS;
}

REF_STATUS ref_sort_unique_int(REF_INT n, REF_INT *original, REF_INT *nunique,
                               REF_INT *unique) {
  REF_INT i, j;

  *nunique = REF_EMPTY;

  RSS(ref_sort_insertion_int(n, original, unique), "sort in unique");

  j = 0;
  for (i = 1; i < n; i++) {
    if (unique[j] != unique[i]) j++;
    if (j != i) unique[j] = unique[i];
  }

  *nunique = j + 1;

  return REF_SUCCESS;
}

REF_STATUS ref_sort_same(REF_INT n, REF_INT *list0, REF_INT *list1,
                         REF_BOOL *same) {
  REF_INT n0, *unique0;
  REF_INT n1, *unique1;
  REF_INT i;
  *same = REF_FALSE;
  ref_malloc(unique0, n, REF_INT);
  ref_malloc(unique1, n, REF_INT);
  RSS(ref_sort_unique_int(n, list0, &n0, unique0), "uniq0");
  RSS(ref_sort_unique_int(n, list1, &n1, unique1), "uniq1");
  if (n0 == n1) {
    *same = REF_TRUE;
    for (i = 0; i < n0; i++) {
      if (unique0[i] != unique1[i]) {
        *same = REF_FALSE;
        break;
      }
    }
  }
  ref_free(unique1);
  ref_free(unique0);
  return REF_SUCCESS;
}

REF_STATUS ref_sort_search_int(REF_INT n, REF_INT *ascending_list,
                               REF_INT target, REF_INT *position) {
  REF_INT lower, upper, mid;

  *position = REF_EMPTY;

  if (n < 1) return REF_NOT_FOUND;

  if (target < ascending_list[0] || target > ascending_list[n - 1])
    return REF_NOT_FOUND;

  lower = 0;
  upper = n - 1;
  mid = n >> 1; /* fast divide by two */

  if (target == ascending_list[lower]) {
    *position = lower;
    return REF_SUCCESS;
  }
  if (target == ascending_list[upper]) {
    *position = upper;
    return REF_SUCCESS;
  }

  while ((lower < mid) && (mid < upper)) {
    if (target >= ascending_list[mid]) {
      if (target == ascending_list[mid]) {
        *position = mid;
        return REF_SUCCESS;
      }
      lower = mid;
    } else {
      upper = mid;
    }
    mid = (lower + upper) >> 1;
  }

  return REF_NOT_FOUND;
}

REF_STATUS ref_sort_search_glob(REF_INT n, REF_GLOB *ascending_list,
                                REF_GLOB target, REF_INT *position) {
  REF_INT lower, upper, mid;

  *position = REF_EMPTY;

  if (n < 1) return REF_NOT_FOUND;

  if (target < ascending_list[0] || target > ascending_list[n - 1])
    return REF_NOT_FOUND;

  lower = 0;
  upper = n - 1;
  mid = n >> 1; /* fast divide by two */

  if (target == ascending_list[lower]) {
    *position = lower;
    return REF_SUCCESS;
  }
  if (target == ascending_list[upper]) {
    *position = upper;
    return REF_SUCCESS;
  }

  while ((lower < mid) && (mid < upper)) {
    if (target >= ascending_list[mid]) {
      if (target == ascending_list[mid]) {
        *position = mid;
        return REF_SUCCESS;
      }
      lower = mid;
    } else {
      upper = mid;
    }
    mid = (lower + upper) >> 1;
  }

  return REF_NOT_FOUND;
}

REF_INT ref_sort_rand_in_range(REF_INT min, REF_INT max) {
  return rand() % (max - min + 1) + min;
}

/* Durstenfeld, R. (July 1964). "Algorithm 235: Random permutation".
 * Communications of the ACM. 7 (7): 420. doi:10.1145/364520.364540
 * Fisher, Ronald A.; Yates, Frank (1948) [1938].
 * Statistical tables for biological, agricultural and medical research
 * (3rd ed.). London: Oliver & Boyd. pp. 26–27. */
REF_STATUS ref_sort_shuffle(REF_INT n, REF_INT *permutation) {
  REF_INT i, j, temp;
  for (i = 0; i < n; i++) permutation[i] = i;
  for (i = 0; i < n - 1; i++) {
    j = rand() % (n - 1 - i + 1) + i; /* i <= j < n */
    j = MAX(i, MIN(j, n - 1));        /* not needed, but don't trust myself */
    temp = permutation[j];
    permutation[j] = permutation[i];
    permutation[i] = temp;
  }
  return REF_SUCCESS;
}
