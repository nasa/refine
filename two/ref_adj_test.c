#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_adj.h"
#include "ref_test.h"

int main( int argc, char *argv[] )
{
  REF_ADJ ref_adj;
  REF_INT degree, item, ref;

  if (argc>1) {printf("%s ignored\n",argv[0]);}

  TFS(ref_adj_free(NULL),"dont free NULL");
  TSS(ref_adj_create(&ref_adj),"create");
  TSS(ref_adj_free(ref_adj),"free");

  /* register and count*/
  TSS(ref_adj_create(&ref_adj),"create");

  TAS(!ref_adj_valid(ref_adj_first(ref_adj,0)),"empty");

  TSS(ref_adj_reg(ref_adj,0,12),"reg");

  item = ref_adj_first(ref_adj,0);
  TES(12,ref_adj_ref(ref_adj,item),"registered ref");

  TSS(ref_adj_free(ref_adj),"free");

  /* iterate */
  TSS(ref_adj_create(&ref_adj),"create");

  degree = 0;
  ref_adj_for(ref_adj,0,item,ref)
    degree++;
  TES(0,degree,"empty degree");

  TSS(ref_adj_reg(ref_adj,0,14),"reg");

  degree = 0;
  ref_adj_for(ref_adj,0,item,ref)
    {
      degree++;
      TES(14,ref,"registered ref");        
    }  
  TES(1,degree,"registered degree");

  TSS(ref_adj_free(ref_adj),"free");

  return 0;
}
