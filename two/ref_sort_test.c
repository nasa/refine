#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_sort.h"

#include "ref_test.h"

int main( void )
{

  {  /* insert sort ordered */
    REF_INT n=4,original[4], sorted[4];
    original[0]=1;
    original[1]=2;
    original[2]=3;
    original[3]=4;
    TSS( ref_sort_insertion( n, original, sorted ), "sort" );
    TEIS( 1, sorted[0], "sorted[0]");
    TEIS( 2, sorted[1], "sorted[1]");
    TEIS( 3, sorted[2], "sorted[2]");
    TEIS( 4, sorted[3], "sorted[3]");
  }

  {  /* insert sort flip */
    REF_INT n=4,original[4], sorted[4];
    original[0]=4;
    original[1]=3;
    original[2]=2;
    original[3]=1;
    TSS( ref_sort_insertion( n, original, sorted ), "sort" );
    TEIS( 1, sorted[0], "sorted[0]");
    TEIS( 2, sorted[1], "sorted[1]");
    TEIS( 3, sorted[2], "sorted[2]");
    TEIS( 4, sorted[3], "sorted[3]");
  }

  {  /* insert sort flip flip */
    REF_INT n=4,original[4], sorted[4];
    original[0]=2;
    original[1]=1;
    original[2]=4;
    original[3]=3;
    TSS( ref_sort_insertion( n, original, sorted ), "sort" );
    TEIS( 1, sorted[0], "sorted[0]");
    TEIS( 2, sorted[1], "sorted[1]");
    TEIS( 3, sorted[2], "sorted[2]");
    TEIS( 4, sorted[3], "sorted[3]");
  }

  {  /* unique */
    REF_INT n=4,original[4], m, unique[4];
    original[0]=2;
    original[1]=1;
    original[2]=2;
    original[3]=3;
    TSS( ref_sort_unique( n, original, &m, unique ), "unique" );
    TEIS( 3, m, "m");
    TEIS( 1, unique[0], "unique[0]");
    TEIS( 2, unique[1], "unique[1]");
    TEIS( 3, unique[2], "unique[2]");
  }

  {  /* search */
    REF_INT n=4,ascending_list[4], position;
    ascending_list[0]=10;
    ascending_list[1]=20;
    ascending_list[2]=30;
    ascending_list[3]=40;

    TSS(ref_sort_search(n,ascending_list,ascending_list[0],&position),"search");
    TEIS( 0, position, "0");
    TSS(ref_sort_search(n,ascending_list,ascending_list[1],&position),"search");
    TEIS( 1, position, "1");
    TSS(ref_sort_search(n,ascending_list,ascending_list[2],&position),"search");
    TEIS( 2, position, "2");
    TSS(ref_sort_search(n,ascending_list,ascending_list[3],&position),"search");
    TEIS( 3, position, "3");

    TEIS(REF_NOT_FOUND, ref_sort_search(n,ascending_list,0,&position),"search");
    TEIS( REF_EMPTY, position, "0");
    TEIS(REF_NOT_FOUND,ref_sort_search(n,ascending_list,15,&position),"search");
    TEIS( REF_EMPTY, position, "15");
    TEIS(REF_NOT_FOUND,ref_sort_search(n,ascending_list,50,&position),"search");
    TEIS( REF_EMPTY, position, "50");
  }

  {  /* search 0 */
    REF_INT n=0, *ascending_list=NULL, position;

    position = 0;
    TEIS(REF_NOT_FOUND,ref_sort_search(n,ascending_list,17,&position),"search");
    TEIS( REF_EMPTY, position, "0");
  }

  {  /* heap sort zero */
    REF_INT n=0,original[1], sorted_index[1];
    original[0]=1;
    sorted_index[0]=2;
    TSS( ref_sort_heap( n, original, sorted_index ), "sort" );
    TEIS( 2, sorted_index[0], "sorted_index[0]");
  }

  {  /* heap sort one */
    REF_INT n=1,original[1], sorted_index[1];
    original[0]=1;
    sorted_index[0]=2;
    TSS( ref_sort_heap( n, original, sorted_index ), "sort" );
    TEIS( 0, sorted_index[0], "sorted_index[0]");
  }

  {  /* heap two order */
    REF_INT n=2,original[2], sorted_index[2];
    original[0]=1;
    original[1]=2;
    TSS( ref_sort_heap( n, original, sorted_index ), "sort" );
    TEIS( 0, sorted_index[0], "sorted_index[0]");
    TEIS( 1, sorted_index[1], "sorted_index[1]");
  }

  {  /* heap two reversed */
    REF_INT n=2,original[2], sorted_index[2];
    original[0]=2;
    original[1]=1;
    TSS( ref_sort_heap( n, original, sorted_index ), "sort" );
    TEIS( 1, sorted_index[0], "sorted_index[0]");
    TEIS( 0, sorted_index[1], "sorted_index[1]");
  }

  {  /* heap three 012 */
    REF_INT n=3,original[3], sorted_index[3];
    original[0]=0;
    original[1]=1;
    original[2]=2;
    TSS( ref_sort_heap( n, original, sorted_index ), "sort" );
    TEIS( 0, sorted_index[0], "sorted_index[0]");
    TEIS( 1, sorted_index[1], "sorted_index[1]");
    TEIS( 2, sorted_index[2], "sorted_index[2]");
  }

  {  /* heap three 012 */
    REF_INT n=3,original[3], sorted_index[3];
    original[0]=0;
    original[1]=1;
    original[2]=2;
    TSS( ref_sort_heap( n, original, sorted_index ), "sort" );
    TEIS( 0, sorted_index[0], "sorted_index[0]");
    TEIS( 1, sorted_index[1], "sorted_index[1]");
    TEIS( 2, sorted_index[2], "sorted_index[2]");
  }

  {  /* heap three 120 */
    REF_INT n=3,original[3], sorted_index[3];
    original[0]=2;
    original[1]=0;
    original[2]=1;
    TSS( ref_sort_heap( n, original, sorted_index ), "sort" );
    TEIS( 1, sorted_index[0], "sorted_index[0]");
    TEIS( 2, sorted_index[1], "sorted_index[1]");
    TEIS( 0, sorted_index[2], "sorted_index[2]");
  }

  return 0;
}
