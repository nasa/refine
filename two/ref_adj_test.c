#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_adj.h"
#include "ref_test.h"

int main( void )
{
  REF_ADJ ref_adj;
  REF_INT degree, item, ref, node;
  REF_INT nitem;

  {
    REIS(REF_NULL, ref_adj_free(NULL),"dont free NULL");
    RSS(ref_adj_create(&ref_adj),"create");
    RSS(ref_adj_free(ref_adj),"free");
  }

  { /* add and count */
    RSS(ref_adj_create(&ref_adj),"create");

    RAS(!ref_adj_valid(ref_adj_first(ref_adj,0)),"empty");

    RSS(ref_adj_add(ref_adj,0,12),"add");

    item = ref_adj_first(ref_adj,0);
    RES(12,ref_adj_safe_ref(ref_adj,item),"added ref");

    RSS(ref_adj_free(ref_adj),"free");
  }

  { /* remove*/
    RSS(ref_adj_create(&ref_adj),"create");

    RSS(ref_adj_add(ref_adj,0,12),"add");
    REIS(REF_INVALID,ref_adj_remove(ref_adj,0,13),"remove missing");
    RSS(ref_adj_remove(ref_adj,0,12),"remove added");

    item = ref_adj_first(ref_adj,0);
    RES(REF_EMPTY,ref_adj_safe_ref(ref_adj,item),"added ref");

    RSS(ref_adj_free(ref_adj),"free");
  }

  { /* iterate */
    RSS(ref_adj_create(&ref_adj),"create");

    degree = 0;
    each_ref_adj_node_item_with_ref(ref_adj,0,item,ref)
      degree++;
    RES(0,degree,"empty degree");

    RSS(ref_adj_add(ref_adj,0,14),"add");

    degree = 0;
    each_ref_adj_node_item_with_ref(ref_adj,0,item,ref)
      {
	degree++;
	RES(14,ref,"check ref");        
      }
    RES(1,degree,"node degree");

    RSS(ref_adj_free(ref_adj),"free");
  }

  { /* empty */
    RSS(ref_adj_create(&ref_adj),"create");

    RAS( ref_adj_empty(ref_adj,0), "starts empty");
    RSS(ref_adj_add(ref_adj,0,14),"add");
    RAS( !ref_adj_empty(ref_adj,0), "not empty anymore");

    RSS(ref_adj_free(ref_adj),"free");
  }

  {  /* negative node */
    RSS(ref_adj_create(&ref_adj),"create");
  
    RES(REF_EMPTY,ref_adj_first(ref_adj,-1),"negative first");
    REIS(REF_INVALID,ref_adj_add(ref_adj,-1,21),"negative add");

    RSS(ref_adj_free(ref_adj),"free");
  }

  { /* reallocate nodes */
    RSS(ref_adj_create(&ref_adj),"create");
    node = ref_adj_nnode( ref_adj );

    RSS(ref_adj_add(ref_adj,node,15),"add and realloc nodes");

    RAS( ref_adj_nnode( ref_adj ) > node, "nodes bigger" );

    item = ref_adj_first(ref_adj,node);
    RES(15,ref_adj_safe_ref(ref_adj,item),"realloc has ref");

    item = ref_adj_first(ref_adj,node+1);
    RES(REF_EMPTY,ref_adj_safe_ref(ref_adj,item),"realloc init empty");

    RSS(ref_adj_free(ref_adj),"free");
  }

  { /* reallocate adj */
    RSS(ref_adj_create(&ref_adj),"create");
    nitem =  ref_adj_nitem( ref_adj );
    for ( item = 0 ; item < nitem+1 ; item++ )
      RSS(ref_adj_add(ref_adj,0,item),"add requiring item realloc");

    RAS( ref_adj_nitem( ref_adj ) > nitem, "item bigger" );

    RSS(ref_adj_free(ref_adj),"free");
  }

  return 0;
}
