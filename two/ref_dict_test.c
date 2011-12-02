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
    TSS(ref_dict_free(ref_dict),"free");
  }

  return 0;
}
