#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_dict.h"
#include "ref_test.h"

int main( void )
{
  REF_DICT ref_dict;

  {
    TFS(ref_dict_free(NULL),"dont free NULL");
    TSS(ref_dict_create(&ref_dict),"create");
    TEIS(0,ref_dict_n(ref_dict),"init zero");
    TSS(ref_dict_free(ref_dict),"free");
  }

  {
    REF_INT key, value;
    TSS(ref_dict_create(&ref_dict),"create");
    key = 2; value = 5;
    TSS(ref_dict_store(ref_dict,key,value),"store");
    TSS(ref_dict_value(ref_dict,key,&value),"retrieve");
    TEIS(5,value,"get value");
    TSS(ref_dict_free(ref_dict),"free");
  }

  return 0;
}
