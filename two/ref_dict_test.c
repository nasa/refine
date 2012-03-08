#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_dict.h"


int main( void )
{
  REF_DICT ref_dict;

  {
    REIS(REF_NULL,ref_dict_free(NULL),"dont free NULL");
    RSS(ref_dict_create(&ref_dict),"create");
    REIS(0,ref_dict_n(ref_dict),"init zero");
    RSS(ref_dict_free(ref_dict),"free");
  }

  { /* missing fails */
    REF_INT key, value;
     key = 2; value = 5;
    RSS(ref_dict_create(&ref_dict),"create");
    REIS(REF_NOT_FOUND,ref_dict_value(ref_dict,key,&value),"missing");
    REIS(5,value,"get value");
    RSS(ref_dict_free(ref_dict),"free");
  }

  { /* store one */
    REF_INT key, value;
    RSS(ref_dict_create(&ref_dict),"create");
    key = 2; value = 5;
    RSS(ref_dict_store(ref_dict,key,value),"store");
    RSS(ref_dict_value(ref_dict,key,&value),"retrieve");
    REIS(5,value,"get value");
    RSS(ref_dict_free(ref_dict),"free");
  }

  { /* store two */
    REF_INT key, value;
    RSS(ref_dict_create(&ref_dict),"create");
    key = 2; value = 5;
    RSS(ref_dict_store(ref_dict,key,value),"store");
    key = 1; value = 3;
    RSS(ref_dict_store(ref_dict,key,value),"store");
    key = 1; RSS(ref_dict_value(ref_dict,key,&value),"retrieve");
    REIS(3,value,"get value");
    key = 2; RSS(ref_dict_value(ref_dict,key,&value),"retrieve");
    REIS(5,value,"get value");
    RSS(ref_dict_free(ref_dict),"free");
  }

  { /* remove */
    REF_INT key, value;
    RSS(ref_dict_create(&ref_dict),"create");
    key = 2; value = 5;
    RSS(ref_dict_store(ref_dict,key,value),"store");
    RSS(ref_dict_remove(ref_dict,key),"remove");
    REIS(0,ref_dict_n(ref_dict),"back to zero");
    REIS(REF_NOT_FOUND,ref_dict_value(ref_dict,key,&value),
	 "should not retrieve");
    RSS(ref_dict_free(ref_dict),"free");
  }

  { /* store lots */
    REF_INT key, value, max;
    RSS(ref_dict_create(&ref_dict),"create");
    max = ref_dict_max(ref_dict);
    for (key=0; key <= max; key++)
      {
	value = 10*key;
	RSS(ref_dict_store(ref_dict,key,value),"store");
      }
    RAS(ref_dict_max(ref_dict)>max, "more?");
    RSS(ref_dict_free(ref_dict),"free");
  }

  { /* store key only once with latest value */
    REF_INT key, value;
    RSS(ref_dict_create(&ref_dict),"create");
    key = 2; value = 5;
    RSS(ref_dict_store(ref_dict,key,value),"store");
    key = 1; value = 3;
    RSS(ref_dict_store(ref_dict,key,value),"store");

    key = 2; value = 7;
    RSS(ref_dict_store(ref_dict,key,value),"store");
    key = 2; RSS(ref_dict_value(ref_dict,key,&value),"retrieve");
    REIS(7,value,"get value");

    REIS(2, ref_dict_n( ref_dict ),"two keys");
    
    RSS(ref_dict_free(ref_dict),"free");
  }

  return 0;
}
