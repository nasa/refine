
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

#include "ref_dict.h"

#include <stdio.h>
#include <stdlib.h>

#include "ref_malloc.h"

REF_STATUS ref_dict_create(REF_DICT *ref_dict_ptr) {
  REF_DICT ref_dict;

  ref_malloc(*ref_dict_ptr, 1, REF_DICT_STRUCT);

  ref_dict = (*ref_dict_ptr);

  ref_dict_n(ref_dict) = 0;
  ref_dict_max(ref_dict) = 10;

  ref_malloc(ref_dict->key, ref_dict_max(ref_dict), REF_INT);
  ref_malloc(ref_dict->value, ref_dict_max(ref_dict), REF_INT);

  return REF_SUCCESS;
}

REF_STATUS ref_dict_free(REF_DICT ref_dict) {
  if (NULL == (void *)ref_dict) return REF_NULL;
  ref_free(ref_dict->value);
  ref_free(ref_dict->key);
  ref_free(ref_dict);
  return REF_SUCCESS;
}

REF_STATUS ref_dict_deep_copy(REF_DICT *ref_dict_ptr, REF_DICT original) {
  REF_DICT ref_dict;
  REF_INT key_index, dict_key, dict_value;

  ref_malloc(*ref_dict_ptr, 1, REF_DICT_STRUCT);

  ref_dict = (*ref_dict_ptr);

  ref_dict_n(ref_dict) = ref_dict_n(original);
  ref_dict_max(ref_dict) = ref_dict_max(original);

  ref_malloc(ref_dict->key, ref_dict_max(ref_dict), REF_INT);
  ref_malloc(ref_dict->value, ref_dict_max(ref_dict), REF_INT);

  each_ref_dict_key_value(original, key_index, dict_key, dict_value) {
    ref_dict_key(ref_dict, key_index) = dict_key;
    ref_dict_keyvalue(ref_dict, key_index) = dict_value;
  }

  return REF_SUCCESS;
}

REF_STATUS ref_dict_store(REF_DICT ref_dict, REF_INT key, REF_INT value) {
  REF_INT i, insert_point;

  if (ref_dict_max(ref_dict) == ref_dict_n(ref_dict)) {
    ref_dict_max(ref_dict) += 1000;

    ref_realloc(ref_dict->key, ref_dict_max(ref_dict), REF_INT);
    ref_realloc(ref_dict->value, ref_dict_max(ref_dict), REF_INT);
  }

  insert_point = 0;
  for (i = ref_dict_n(ref_dict) - 1; i >= 0; i--) {
    if (ref_dict->key[i] == key) {
      ref_dict->value[i] = value;
      return REF_SUCCESS;
    }
    if (ref_dict->key[i] < key) {
      insert_point = i + 1;
      break;
    }
  }
  /* shift to open up insert_point */
  for (i = ref_dict_n(ref_dict); i > insert_point; i--)
    ref_dict->key[i] = ref_dict->key[i - 1];
  for (i = ref_dict_n(ref_dict); i > insert_point; i--)
    ref_dict->value[i] = ref_dict->value[i - 1];

  /* fill insert_point */
  ref_dict_n(ref_dict)++;
  ref_dict->key[insert_point] = key;
  ref_dict->value[insert_point] = value;

  return REF_SUCCESS;
}

REF_STATUS ref_dict_location(REF_DICT ref_dict, REF_INT key,
                             REF_INT *location) {
  REF_INT i;

  *location = REF_EMPTY;

  for (i = 0; i < ref_dict_n(ref_dict); i++)
    if (key == ref_dict->key[i]) {
      *location = i;
      return REF_SUCCESS;
    }

  return REF_NOT_FOUND;
}

REF_STATUS ref_dict_remove(REF_DICT ref_dict, REF_INT key) {
  REF_INT i, location;

  RAISE(ref_dict_location(ref_dict, key, &location));

  ref_dict_n(ref_dict)--;

  for (i = location; i < ref_dict_n(ref_dict); i++) {
    ref_dict->key[i] = ref_dict->key[i + 1];
  }
  for (i = location; i < ref_dict_n(ref_dict); i++) {
    ref_dict->value[i] = ref_dict->value[i + 1];
  }

  return REF_SUCCESS;
}

REF_STATUS ref_dict_value(REF_DICT ref_dict, REF_INT key, REF_INT *value) {
  REF_INT location;

  RAISE(ref_dict_location(ref_dict, key, &location));

  *value = ref_dict->value[location];

  return REF_SUCCESS;
}

REF_BOOL ref_dict_has_key(REF_DICT ref_dict, REF_INT key) {
  REF_INT i;

  for (i = 0; i < ref_dict_n(ref_dict); i++)
    if (key == ref_dict->key[i]) {
      return REF_TRUE;
    }

  return REF_FALSE;
}

REF_BOOL ref_dict_has_value(REF_DICT ref_dict, REF_INT value) {
  REF_INT i;

  for (i = 0; i < ref_dict_n(ref_dict); i++) {
    if (value == ref_dict->value[i]) {
      return REF_TRUE;
    }
  }

  return REF_FALSE;
}

REF_STATUS ref_dict_inspect(REF_DICT ref_dict) {
  REF_INT i;

  printf("ref_dict = %p\n", (void *)ref_dict);
  printf(" n = %d, max = %d\n", ref_dict_n(ref_dict), ref_dict_max(ref_dict));

  each_ref_dict_key_index(ref_dict, i) {
    printf(" %d [%d] = %d\n", ref_dict->key[i], i, ref_dict->value[i]);
  }

  return REF_SUCCESS;
}

REF_STATUS ref_dict_inspect_keys(REF_DICT ref_dict) {
  REF_INT i;

  each_ref_dict_key_index(ref_dict, i) {
    printf(" %d: %d\n", i, ref_dict->key[i]);
  }

  return REF_SUCCESS;
}
