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

  RSS( ref_mpi_start( argc, argv ), "start" );

  {
    REF_INT last;
    TFS(ref_list_free(NULL),"dont free NULL");
    RSS(ref_list_create(&ref_list),"create");
    REIS(0,ref_list_n(ref_list),"init zero");
    TFS(ref_list_remove(ref_list,&last),"rm");
    REIS(REF_EMPTY,last,"remove empty");
    RSS(ref_list_free(ref_list),"free");
  }

  { /* store one */
    REF_INT item;
    RSS(ref_list_create(&ref_list),"create");
    item = 27;
    RSS(ref_list_add(ref_list,item),"add");
    REIS(1,ref_list_n(ref_list),"has one");
    RSS(ref_list_free(ref_list),"free");
  }

  { /* remove */
    REF_INT item, last;
    RSS(ref_list_create(&ref_list),"create");
    item = 27;
    RSS(ref_list_add(ref_list,item),"add");
    RSS(ref_list_remove(ref_list,&last),"rm");
    REIS(0,ref_list_n(ref_list),"has none");

    RSS(ref_list_free(ref_list),"free");
  }

  { /* store lots */
    REF_INT item, max;
    RSS(ref_list_create(&ref_list),"create");
    max = ref_list_max(ref_list);
    for (item=0; item <= max; item++)
      {
	RSS(ref_list_add(ref_list,item),"store");
      }
    TAS(ref_list_max(ref_list)>max, "more?");
    RSS(ref_list_free(ref_list),"free");
  }

  { /* shift */
    REF_INT last;
    RSS(ref_list_create(&ref_list),"create");
    RSS(ref_list_add(ref_list,20),"store");
    RSS(ref_list_add(ref_list,10),"store");
    RSS(ref_list_shift(ref_list,15,27),"shift");

    RSS(ref_list_remove(ref_list,&last),"rm");
    REIS(10,last,"has none");

    RSS(ref_list_remove(ref_list,&last),"rm");
    REIS(47,last,"has none");

    RSS(ref_list_free(ref_list),"free");
  }

  { /* sort */
    REF_INT last;
    RSS(ref_list_create(&ref_list),"create");
    RSS(ref_list_add(ref_list,20),"store");
    RSS(ref_list_add(ref_list,10),"store");
    RSS(ref_list_sort(ref_list),"sort");

    RSS(ref_list_remove(ref_list,&last),"rm");
    REIS(20,last,"has none");
    RSS(ref_list_remove(ref_list,&last),"rm");
    REIS(10,last,"has none");

    RSS(ref_list_free(ref_list),"free");
  }

  { /* erase */
    RSS(ref_list_create(&ref_list),"create");
    RSS(ref_list_add(ref_list,20),"store");
    RSS(ref_list_add(ref_list,10),"store");
    RSS(ref_list_sort(ref_list),"sort");

    RSS(ref_list_erase(ref_list),"rm -rf");

    REIS(0,ref_list_n(ref_list),"has none");

    RSS(ref_list_free(ref_list),"free");
  }

  { /* allgather */
    RSS(ref_list_create(&ref_list),"create");

    RSS(ref_list_add(ref_list,ref_mpi_id),"store");

    RSS(ref_list_allgather(ref_list),"gather");

    REIS(ref_mpi_n,ref_list_n(ref_list),"one from each");

    RSS(ref_list_free(ref_list),"free");
  }

  RSS( ref_mpi_stop( ), "stop" );

  return 0;
}
