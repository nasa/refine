
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

#include "ref_adj.h"

#include <stdio.h>
#include <stdlib.h>

#include "ref_malloc.h"

REF_STATUS ref_adj_create(REF_ADJ *ref_adj_ptr) {
  REF_ADJ ref_adj;
  REF_INT i;

  ref_malloc(*ref_adj_ptr, 1, REF_ADJ_STRUCT);
  ref_adj = (*ref_adj_ptr);

  ref_adj_nnode(ref_adj) = 10;
  ref_adj_nitem(ref_adj) = 20;

  ref_malloc(ref_adj->first, ref_adj_nnode(ref_adj), REF_INT);
  for (i = 0; i < ref_adj_nnode(ref_adj); i++) ref_adj->first[i] = REF_EMPTY;

  ref_malloc(ref_adj->item, ref_adj_nitem(ref_adj), REF_ADJ_ITEM_STRUCT);
  for (i = 0; i < ref_adj_nitem(ref_adj); i++) {
    ref_adj->item[i].ref = REF_EMPTY;
    ref_adj->item[i].next = i + 1;
  }
  ref_adj->item[ref_adj_nitem(ref_adj) - 1].next = REF_EMPTY;
  ref_adj->blank = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_adj_free(REF_ADJ ref_adj) {
  if (NULL == (void *)ref_adj) return REF_NULL;
  ref_free(ref_adj->first);
  ref_free(ref_adj->item);
  ref_free(ref_adj);
  return REF_SUCCESS;
}

REF_STATUS ref_adj_deep_copy(REF_ADJ *ref_adj_ptr, REF_ADJ original) {
  REF_ADJ ref_adj;
  REF_INT i;

  ref_malloc(*ref_adj_ptr, 1, REF_ADJ_STRUCT);
  ref_adj = (*ref_adj_ptr);

  ref_adj_nnode(ref_adj) = ref_adj_nnode(original);
  ref_adj_nitem(ref_adj) = ref_adj_nitem(original);

  ref_malloc(ref_adj->first, ref_adj_nnode(ref_adj), REF_INT);
  for (i = 0; i < ref_adj_nnode(ref_adj); i++)
    ref_adj->first[i] = original->first[i];

  ref_malloc(ref_adj->item, ref_adj_nitem(ref_adj), REF_ADJ_ITEM_STRUCT);
  for (i = 0; i < ref_adj_nitem(ref_adj); i++) {
    ref_adj->item[i].ref = original->item[i].ref;
    ref_adj->item[i].next = original->item[i].next;
  }
  ref_adj->blank = original->blank;

  return REF_SUCCESS;
}

REF_STATUS ref_adj_inspect(REF_ADJ ref_adj) {
  REF_INT node, item;
  printf("ref_adj = %p\n", (void *)ref_adj);
  printf(" blank = %d\n", ref_adj->blank);
  for (node = 0; node < ref_adj_nnode(ref_adj); node++)
    printf(" first[%d] = %d\n", node, ref_adj->first[node]);
  for (item = 0; item < ref_adj_nitem(ref_adj); item++)
    printf(" item[%d] = %d : %d\n", item, ref_adj->item[item].next,
           ref_adj->item[item].ref);

  return REF_SUCCESS;
}

REF_STATUS ref_adj_node_inspect(REF_ADJ ref_adj, REF_INT node) {
  REF_INT item, ref;
  printf(" %d :", node);
  each_ref_adj_node_item_with_ref(ref_adj, node, item, ref) {
    printf(" %d", ref);
  }
  printf("\n");

  return REF_SUCCESS;
}

REF_STATUS ref_adj_add(REF_ADJ ref_adj, REF_INT node, REF_INT reference) {
  REF_INT item;
  REF_INT orig, chunk, i;
  REF_INT max_limit = REF_INT_MAX;

  if (node < 0) return REF_INVALID;

  if (node >= ref_adj_nnode(ref_adj)) {
    orig = ref_adj_nnode(ref_adj);
    chunk = 100 + MAX(0, node - orig);
    chunk = MAX(chunk, (REF_INT)(0.5 * (REF_DBL)orig));
    /* try to keep under 32-bit limit */
    chunk = MIN(chunk, max_limit - orig);
    ref_adj_nnode(ref_adj) = orig + chunk;
    ref_realloc(ref_adj->first, ref_adj_nnode(ref_adj), REF_INT);
    for (i = orig; i < ref_adj_nnode(ref_adj); i++)
      ref_adj->first[i] = REF_EMPTY;
  }

  if (REF_EMPTY == ref_adj_blank(ref_adj)) {
    RAS(ref_adj_nitem(ref_adj) != max_limit,
        "the number of ref_adj items is too large for int, cannot grow");

    orig = ref_adj_nitem(ref_adj);
    chunk = MAX(100, (REF_INT)(0.5 * (REF_DBL)orig));
    /* try to keep under 32-bit limit */
    chunk = MIN(chunk, max_limit - orig);
    ref_adj_nitem(ref_adj) = orig + chunk;
    ref_realloc(ref_adj->item, ref_adj_nitem(ref_adj), REF_ADJ_ITEM_STRUCT);
    for (i = orig; i < ref_adj_nitem(ref_adj); i++) {
      ref_adj->item[i].ref = REF_EMPTY;
      ref_adj->item[i].next = i + 1;
    }
    ref_adj->item[ref_adj_nitem(ref_adj) - 1].next = REF_EMPTY;
    ref_adj->blank = orig;
  }

  item = ref_adj_blank(ref_adj);
  ref_adj_blank(ref_adj) = ref_adj_item_next(ref_adj, item);

  ref_adj_item_ref(ref_adj, item) = reference;
  ref_adj_item_next(ref_adj, item) = ref_adj_first(ref_adj, node);

  ref_adj->first[node] = item;

  return REF_SUCCESS;
}

REF_STATUS ref_adj_remove(REF_ADJ ref_adj, REF_INT node, REF_INT reference) {
  REF_INT item, ref;
  REF_INT target, parent;

  item = ref_adj_first(ref_adj, node);

  if (!ref_adj_valid(item)) return REF_INVALID;

  if (reference == ref_adj_item_ref(ref_adj, item)) {
    ref_adj->first[node] = ref_adj_item_next(ref_adj, item);
    ref_adj_item_next(ref_adj, item) = ref_adj_blank(ref_adj);
    ref_adj_blank(ref_adj) = item;
    ref_adj_item_ref(ref_adj, item) = REF_EMPTY;
    return REF_SUCCESS;
  }

  target = REF_EMPTY;
  parent = REF_EMPTY;
  each_ref_adj_node_item_with_ref(ref_adj, node, item, ref) {
    if (ref == reference) {
      target = item;
      break;
    } else {
      parent = item;
    }
  }

  if (REF_EMPTY == target) return REF_INVALID;

  if (REF_EMPTY == parent) RSS(REF_FAILURE, "parent empty");

  ref_adj_item_next(ref_adj, parent) = ref_adj_item_next(ref_adj, item);
  ref_adj_item_next(ref_adj, item) = ref_adj_blank(ref_adj);
  ref_adj_blank(ref_adj) = item;
  ref_adj_item_ref(ref_adj, item) = REF_EMPTY;

  return REF_SUCCESS;
}

REF_STATUS ref_adj_add_uniquely(REF_ADJ ref_adj, REF_INT node,
                                REF_INT reference) {
  REF_INT item, ref;

  each_ref_adj_node_item_with_ref(ref_adj, node, item, ref) {
    if (ref == reference) {
      return REF_SUCCESS;
    }
  }

  return ref_adj_add(ref_adj, node, reference);
}

REF_STATUS ref_adj_degree(REF_ADJ ref_adj, REF_INT node, REF_INT *degree) {
  REF_INT item;
  *degree = 0;

  for (item = ref_adj_first(ref_adj, node); ref_adj_valid(item);
       (item) = ref_adj_item_next(ref_adj, item))
    (*degree)++;

  return REF_SUCCESS;
}

REF_STATUS ref_adj_min_degree_node(REF_ADJ ref_adj, REF_INT *min_degree,
                                   REF_INT *min_degree_node) {
  REF_INT node, degree;
  *min_degree = REF_EMPTY;
  *min_degree_node = REF_EMPTY;

  for (node = 0; node < ref_adj_nnode(ref_adj); node++) {
    RSS(ref_adj_degree(ref_adj, node, &degree), "deg");
    if (degree > 0) {
      if (REF_EMPTY == (*min_degree_node) || degree < (*min_degree)) {
        *min_degree_node = node;
        *min_degree = degree;
      }
    }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_adj_tec_fill(REF_ADJ ref_adj, const char *filename) {
  REF_INT node, item, nadj;

  FILE *file;

  file = fopen(filename, "w");
  if (NULL == (void *)file) printf("unable to open %s\n", filename);
  RNS(file, "unable to open file");

  fprintf(file, "title=\"tecplot refine adj fill\"\n");
  fprintf(file, "variables = \"item\" \"node\"\n");

  nadj = 0;
  for (node = 0; node < ref_adj_nnode(ref_adj); node++) {
    each_ref_adj_node_item(ref_adj, node, item) { nadj++; }
  }

  fprintf(file, "zone t=\"fill\", i=%d, datapacking=%s\n", nadj, "point");

  for (node = 0; node < ref_adj_nnode(ref_adj); node++) {
    each_ref_adj_node_item(ref_adj, node, item) {
      fprintf(file, " %d %d\n", item, node);
    }
  }

  fclose(file);

  return REF_SUCCESS;
}
