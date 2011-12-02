
#ifndef REF_DICT_H
#define REF_DICT_H

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_DICT_STRUCT REF_DICT_STRUCT;
typedef REF_DICT_STRUCT * REF_DICT;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_DICT_STRUCT {
  REF_INT n, max;
  REF_INT *key;
  REF_INT *value;
};

REF_STATUS ref_dict_create( REF_DICT *ref_dict );
REF_STATUS ref_dict_free( REF_DICT ref_dict );

#define ref_dict_n( ref_dict ) ((ref_dict)->n)
#define ref_dict_max( ref_dict ) ((ref_dict)->max)

REF_STATUS ref_dict_inspect( REF_DICT ref_dict );

END_C_DECLORATION

#endif /* REF_DICT_H */
