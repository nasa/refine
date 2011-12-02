
#include <stdlib.h>
#include <stdio.h>

#include "ref_dict.h"

REF_STATUS ref_dict_create( REF_DICT *ref_dict_ptr )
{
  REF_DICT ref_dict;

  (*ref_dict_ptr) = NULL;
  (*ref_dict_ptr) = (REF_DICT)malloc( sizeof(REF_DICT_STRUCT) );
  RNS(*ref_dict_ptr,"malloc ref_dict NULL");
  ref_dict = (*ref_dict_ptr);

  ref_dict_n(ref_dict) = 0;
  ref_dict_max(ref_dict) = 10;

  ref_dict->key = (REF_INT *)malloc( ref_dict_max(ref_dict) * 
				     sizeof(REF_INT) );
  RNS(ref_dict->key,"malloc ref_dict->key NULL");

  ref_dict->value = (REF_INT *)malloc( ref_dict_max(ref_dict) * 
				     sizeof(REF_INT) );
  RNS(ref_dict->value,"malloc ref_dict->value NULL");

  return REF_SUCCESS;
}

REF_STATUS ref_dict_free( REF_DICT ref_dict )
{
  if ( NULL == (void *)ref_dict ) return REF_NULL;
  ref_cond_free( ref_dict->value );
  ref_cond_free( ref_dict->key );
  ref_cond_free( ref_dict );
  return REF_SUCCESS;
}

REF_STATUS ref_dict_store( REF_DICT ref_dict, REF_INT key, REF_INT value )
{
  ref_dict->key[0] = key;
  ref_dict->value[0] = value;
  return REF_SUCCESS;
}

REF_STATUS ref_dict_value( REF_DICT ref_dict, REF_INT key, REF_INT *value )
{
  if ( key == ref_dict->key[0] )
    *value = ref_dict->value[0];
  return REF_SUCCESS;
}

REF_STATUS ref_dict_inspect( REF_DICT ref_dict )
{
  REF_INT i;
  printf("ref_dict = %p\n",(void *)ref_dict);
  printf(" n = %d, max = %d\n",ref_dict_n(ref_dict),ref_dict_max(ref_dict));
  for ( i = 0 ; i < ref_dict_n( ref_dict ) ; i++ )
    printf(" %d [%d] = %d\n",ref_dict->key[i],i,ref_dict->value[i]);

  return REF_SUCCESS;
}

