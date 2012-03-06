#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_list.h"
#include "ref_test.h"

int main( void )
{
  REF_LIST ref_list;

  {
    REF_INT last;
    TFS(ref_list_free(NULL),"dont free NULL");
    TSS(ref_list_create(&ref_list),"create");
    TEIS(0,ref_list_n(ref_list),"init zero");
    TFS(ref_list_remove(ref_list,&last),"rm");
    TEIS(REF_EMPTY,last,"remove empty");
    TSS(ref_list_free(ref_list),"free");
  }

  { /* store one */
    REF_INT item;
    TSS(ref_list_create(&ref_list),"create");
    item = 27;
    TSS(ref_list_add(ref_list,item),"add");
    TEIS(1,ref_list_n(ref_list),"has one");
    TSS(ref_list_free(ref_list),"free");
  }

  { /* remove */
    REF_INT item, last;
    TSS(ref_list_create(&ref_list),"create");
    item = 27;
    TSS(ref_list_add(ref_list,item),"add");
    TSS(ref_list_remove(ref_list,&last),"rm");
    TEIS(0,ref_list_n(ref_list),"has none");

    TSS(ref_list_free(ref_list),"free");
  }

  { /* store lots */
    REF_INT item, max;
    TSS(ref_list_create(&ref_list),"create");
    max = ref_list_max(ref_list);
    for (item=0; item <= max; item++)
      {
	TSS(ref_list_add(ref_list,item),"store");
      }
    TAS(ref_list_max(ref_list)>max, "more?");
    TSS(ref_list_free(ref_list),"free");
  }

  { /* shift */
    REF_INT last;
    TSS(ref_list_create(&ref_list),"create");
    TSS(ref_list_add(ref_list,20),"store");
    TSS(ref_list_add(ref_list,10),"store");
    TSS(ref_list_shift(ref_list,15,27),"shift");

    TSS(ref_list_remove(ref_list,&last),"rm");
    TEIS(10,last,"has none");

    TSS(ref_list_remove(ref_list,&last),"rm");
    TEIS(47,last,"has none");

    TSS(ref_list_free(ref_list),"free");
  }

  return 0;
}
