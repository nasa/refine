#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_adj.h"
#include "ref_test.h"

int main( int argc, char *argv[] )
{
  REF_ADJ ref_adj;

  if (argc>1) {printf("%s ignored\n",argv[0]);}

  TFS(ref_adj_free(NULL),"dont free NULL");

  TSS(ref_adj_create(&ref_adj),"create");

  TSS(ref_adj_free(ref_adj),"free");

  return 0;
}
