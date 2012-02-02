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

  return 0;
}
