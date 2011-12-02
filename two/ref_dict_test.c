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

  { /* missing fails */
    REF_INT key, value;
     key = 2; value = 5;
    TSS(ref_dict_create(&ref_dict),"create");
    TFS(ref_dict_value(ref_dict,key,&value),"missing");
    TEIS(5,value,"get value");
    TSS(ref_dict_free(ref_dict),"free");
  }

  { /* store one */
    REF_INT key, value;
    TSS(ref_dict_create(&ref_dict),"create");
    key = 2; value = 5;
    TSS(ref_dict_store(ref_dict,key,value),"store");
    TSS(ref_dict_value(ref_dict,key,&value),"retrieve");
    TEIS(5,value,"get value");
    TSS(ref_dict_free(ref_dict),"free");
  }

  { /* store two */
    REF_INT key, value;
    TSS(ref_dict_create(&ref_dict),"create");
    key = 2; value = 5;
    TSS(ref_dict_store(ref_dict,key,value),"store");
    key = 1; value = 3;
    TSS(ref_dict_store(ref_dict,key,value),"store");
    key = 1; TSS(ref_dict_value(ref_dict,key,&value),"retrieve");
    TEIS(3,value,"get value");
    key = 2; TSS(ref_dict_value(ref_dict,key,&value),"retrieve");
    TEIS(5,value,"get value");
    TSS(ref_dict_free(ref_dict),"free");
  }

  { /* remove */
    REF_INT key, value;
    TSS(ref_dict_create(&ref_dict),"create");
    key = 2; value = 5;
    TSS(ref_dict_store(ref_dict,key,value),"store");
    TSS(ref_dict_remove(ref_dict,key),"remove");
    TEIS(0,ref_dict_n(ref_dict),"back to zero");
    TFS(ref_dict_value(ref_dict,key,&value),"should not retrieve");
    TSS(ref_dict_free(ref_dict),"free");
  }

  SKIP_TEST("realloc"){};

  return 0;
}
