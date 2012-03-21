#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_sort.h"



int main( void )
{

  {  /* insert sort ordered */
    REF_INT n=4,original[4], sorted[4];
    original[0]=1;
    original[1]=2;
    original[2]=3;
    original[3]=4;
    RSS( ref_sort_insertion_int( n, original, sorted ), "sort" );
    REIS( 1, sorted[0], "sorted[0]");
    REIS( 2, sorted[1], "sorted[1]");
    REIS( 3, sorted[2], "sorted[2]");
    REIS( 4, sorted[3], "sorted[3]");
  }

  {  /* insert sort flip */
    REF_INT n=4,original[4], sorted[4];
    original[0]=4;
    original[1]=3;
    original[2]=2;
    original[3]=1;
    RSS( ref_sort_insertion_int( n, original, sorted ), "sort" );
    REIS( 1, sorted[0], "sorted[0]");
    REIS( 2, sorted[1], "sorted[1]");
    REIS( 3, sorted[2], "sorted[2]");
    REIS( 4, sorted[3], "sorted[3]");
  }

  {  /* insert sort flip flip */
    REF_INT n=4,original[4], sorted[4];
    original[0]=2;
    original[1]=1;
    original[2]=4;
    original[3]=3;
    RSS( ref_sort_insertion_int( n, original, sorted ), "sort" );
    REIS( 1, sorted[0], "sorted[0]");
    REIS( 2, sorted[1], "sorted[1]");
    REIS( 3, sorted[2], "sorted[2]");
    REIS( 4, sorted[3], "sorted[3]");
  }

  {  /* unique */
    REF_INT n=4,original[4], m, unique[4];
    original[0]=2;
    original[1]=1;
    original[2]=2;
    original[3]=3;
    RSS( ref_sort_unique_int( n, original, &m, unique ), "unique" );
    REIS( 3, m, "m");
    REIS( 1, unique[0], "unique[0]");
    REIS( 2, unique[1], "unique[1]");
    REIS( 3, unique[2], "unique[2]");
  }

  {  /* search */
    REF_INT n=4,ascending_list[4], position;
    ascending_list[0]=10;
    ascending_list[1]=20;
    ascending_list[2]=30;
    ascending_list[3]=40;

    RSS(ref_sort_search(n,ascending_list,ascending_list[0],&position),"search");
    REIS( 0, position, "0");
    RSS(ref_sort_search(n,ascending_list,ascending_list[1],&position),"search");
    REIS( 1, position, "1");
    RSS(ref_sort_search(n,ascending_list,ascending_list[2],&position),"search");
    REIS( 2, position, "2");
    RSS(ref_sort_search(n,ascending_list,ascending_list[3],&position),"search");
    REIS( 3, position, "3");

    REIS(REF_NOT_FOUND, ref_sort_search(n,ascending_list,0,&position),"search");
    REIS( REF_EMPTY, position, "0");
    REIS(REF_NOT_FOUND,ref_sort_search(n,ascending_list,15,&position),"search");
    REIS( REF_EMPTY, position, "15");
    REIS(REF_NOT_FOUND,ref_sort_search(n,ascending_list,50,&position),"search");
    REIS( REF_EMPTY, position, "50");
  }

  {  /* search 0 */
    REF_INT n=0, *ascending_list=NULL, position;

    position = 0;
    REIS(REF_NOT_FOUND,ref_sort_search(n,ascending_list,17,&position),"search");
    REIS( REF_EMPTY, position, "0");
  }

  {  /* heap sort zero */
    REF_INT n=0,original[1], sorted_index[1];
    original[0]=1;
    sorted_index[0]=2;
    RSS( ref_sort_heap_int( n, original, sorted_index ), "sort" );
    REIS( 2, sorted_index[0], "sorted_index[0]");
  }

  {  /* heap sort one */
    REF_INT n=1,original[1], sorted_index[1];
    original[0]=1;
    sorted_index[0]=2;
    RSS( ref_sort_heap_int( n, original, sorted_index ), "sort" );
    REIS( 0, sorted_index[0], "sorted_index[0]");
  }

  {  /* heap two order */
    REF_INT n=2,original[2], sorted_index[2];
    original[0]=1;
    original[1]=2;
    RSS( ref_sort_heap_int( n, original, sorted_index ), "sort" );
    REIS( 0, sorted_index[0], "sorted_index[0]");
    REIS( 1, sorted_index[1], "sorted_index[1]");
  }

  {  /* heap two reversed */
    REF_INT n=2,original[2], sorted_index[2];
    original[0]=2;
    original[1]=1;
    RSS( ref_sort_heap_int( n, original, sorted_index ), "sort" );
    REIS( 1, sorted_index[0], "sorted_index[0]");
    REIS( 0, sorted_index[1], "sorted_index[1]");
  }

  {  /* heap three 012 */
    REF_INT n=3,original[3], sorted_index[3];
    original[0]=0;
    original[1]=1;
    original[2]=2;
    RSS( ref_sort_heap_int( n, original, sorted_index ), "sort" );
    REIS( 0, sorted_index[0], "sorted_index[0]");
    REIS( 1, sorted_index[1], "sorted_index[1]");
    REIS( 2, sorted_index[2], "sorted_index[2]");
  }

  {  /* heap three 012 */
    REF_INT n=3,original[3], sorted_index[3];
    original[0]=0;
    original[1]=1;
    original[2]=2;
    RSS( ref_sort_heap_int( n, original, sorted_index ), "sort" );
    REIS( 0, sorted_index[0], "sorted_index[0]");
    REIS( 1, sorted_index[1], "sorted_index[1]");
    REIS( 2, sorted_index[2], "sorted_index[2]");
  }

  {  /* heap three 120 */
    REF_INT n=3,original[3], sorted_index[3];
    original[0]=2;
    original[1]=0;
    original[2]=1;
    RSS( ref_sort_heap_int( n, original, sorted_index ), "sort" );
    REIS( 1, sorted_index[0], "sorted_index[0]");
    REIS( 2, sorted_index[1], "sorted_index[1]");
    REIS( 0, sorted_index[2], "sorted_index[2]");
  }

  return 0;
}
