#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_list.h"
#include "ref_sort.h"
#include "ref_mpi.h"

#include "ref_test.h"

int main( int argc, char *argv[] )
{
  REF_LIST ref_list;

  TSS( ref_mpi_start( argc, argv ), "start" );

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

  { /* sort */
    REF_INT last;
    TSS(ref_list_create(&ref_list),"create");
    TSS(ref_list_add(ref_list,20),"store");
    TSS(ref_list_add(ref_list,10),"store");
    TSS(ref_list_sort(ref_list),"sort");

    TSS(ref_list_remove(ref_list,&last),"rm");
    TEIS(20,last,"has none");
    TSS(ref_list_remove(ref_list,&last),"rm");
    TEIS(10,last,"has none");

    TSS(ref_list_free(ref_list),"free");
  }

  { /* erase */
    TSS(ref_list_create(&ref_list),"create");
    TSS(ref_list_add(ref_list,20),"store");
    TSS(ref_list_add(ref_list,10),"store");
    TSS(ref_list_sort(ref_list),"sort");

    TSS(ref_list_erase(ref_list),"rm -rf");

    TEIS(0,ref_list_n(ref_list),"has none");

    TSS(ref_list_free(ref_list),"free");
  }

  { /* allgather */
    TSS(ref_list_create(&ref_list),"create");

    TSS(ref_list_add(ref_list,ref_mpi_id),"store");

    TSS(ref_list_allgather(ref_list),"gather");

    TEIS(ref_mpi_n,ref_list_n(ref_list),"one from each");

    TSS(ref_list_free(ref_list),"free");
  }

  TSS( ref_mpi_stop( ), "stop" );

  return 0;
}
