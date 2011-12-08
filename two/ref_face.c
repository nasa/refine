
#include <stdlib.h>
#include <stdio.h>

#include "ref_face.h"

REF_STATUS ref_face_create( REF_FACE *ref_face_ptr, REF_GRID ref_grid )
{
  REF_FACE ref_face;
  REF_INT group, cell, cell_face;
  REF_INT node;
  REF_INT nodes[4];
  REF_CELL ref_cell;

  (*ref_face_ptr) = NULL;
  (*ref_face_ptr) = (REF_FACE)malloc( sizeof(REF_FACE_STRUCT) );
  RNS(*ref_face_ptr,"malloc ref_face NULL");

  ref_face = *ref_face_ptr;

  ref_face_n(ref_face) = 0;
  ref_face_max(ref_face) = 10;

  ref_face->f2n = (REF_INT *)malloc( 4 * ref_face_max(ref_face) * 
				     sizeof(REF_INT) );
  RNS(ref_face->f2n,"malloc ref_face->f2n NULL");

  RSS( ref_adj_create( &(ref_face_adj( ref_face )) ), "create adj" );

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    each_ref_cell_valid_cell( ref_cell, cell )
      each_ref_cell_cell_face( ref_cell, cell_face )
        {
	  for(node=0;node<4;node++)
	    nodes[node]=ref_cell_f2n(ref_cell,node,cell,cell_face);
	  RSS( ref_face_add_uniquely( ref_face, nodes ), "add face");
	}

  return REF_SUCCESS;
}

REF_STATUS ref_face_free( REF_FACE ref_face )
{
  if ( NULL == (void *)ref_face ) return REF_NULL;

  RSS( ref_adj_free( ref_face_adj( ref_face ) ), "free adj" );
  free( ref_face->f2n );

  ref_cond_free( ref_face );

  return REF_SUCCESS;
}

REF_STATUS ref_face_inspect( REF_FACE ref_face )
{
  REF_INT face,node,orig[4],sort[4];

  printf("ref_face = %p\n",(void *)ref_face);
  printf(" n = %d\n",ref_face_n(ref_face));
  printf(" max = %d\n",ref_face_max(ref_face));
  for (face=0;face<ref_face_n(ref_face);face++)
    printf("face %4d : %4d %4d %4d %4d\n",face,
	   ref_face_f2n(ref_face,0,face),
	   ref_face_f2n(ref_face,1,face),
	   ref_face_f2n(ref_face,2,face),
	   ref_face_f2n(ref_face,3,face) );
  for (face=0;face<ref_face_n(ref_face);face++)
    {
      for(node=0;node<4;node++)
	orig[node]=ref_face_f2n(ref_face,node,face);
      RSS( ref_insertion_sort( 4, orig, sort ), "sort" );
      printf("sort %4d : %4d %4d %4d %4d\n",face,
	     sort[0],sort[1],sort[2],sort[3]);
    }
  return REF_SUCCESS;
}

REF_STATUS ref_insertion_sort( REF_INT n, REF_INT *original, REF_INT *sorted )
{
  REF_INT i, j, smallest, temp;

  for(i=0;i<n;i++)
    sorted[i] = original[i];

  for(i=0;i<n;i++)
    {
      smallest = i;
      for(j=i+1;j<n;j++)
	{
	  if ( sorted[j] < sorted[smallest] )
	    smallest = j;
	}
      temp = sorted[i];
      sorted[i] = sorted[smallest];
      sorted[smallest] = temp;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_face_with( REF_FACE ref_face, REF_INT *nodes, REF_INT *face )
{
  REF_INT item, ref, node;
  REF_INT target[4], canidate[4], orig[4];

  printf("\n");

  (*face) = REF_EMPTY;

  RSS( ref_insertion_sort( 4, nodes, target ), "sort" );

  printf("target   : %3d %3d %3d %3d\n",
	 target[0],target[1],target[2],target[3]);

  each_ref_adj_node_item_with_ref( ref_face_adj(ref_face), nodes[0], item, ref)
    {
      for(node=0;node<4;node++)
	orig[node]=ref_face_f2n(ref_face,node,ref);
      RSS( ref_insertion_sort( 4, orig, canidate ), "sort" );
  printf("canidate : %3d %3d %3d %3d\n",
	 canidate[0],canidate[1],canidate[2],canidate[3]);
      if ( target[0] == canidate[0] &&
	   target[1] == canidate[1] &&
	   target[2] == canidate[2] &&
	   target[3] == canidate[3] )
	{
	  printf("found\n");
	  (*face) = ref;
	  return REF_SUCCESS;
	}
    }

  return REF_NOT_FOUND;
}

REF_STATUS ref_face_add_uniquely( REF_FACE ref_face, REF_INT *nodes )
{
  REF_INT face, node;
  REF_STATUS status;

  status = ref_face_with( ref_face, nodes, &face );
  if ( REF_SUCCESS == status ) return REF_SUCCESS;
  if ( REF_NOT_FOUND != status ) RSS( status, "looking for face");

  if ( ref_face_n(ref_face) >= ref_face_max(ref_face) )
    {
      REF_INT chunk;
      chunk = 1000;
      ref_face_max(ref_face) += chunk;
      ref_face->f2n = (REF_INT *)realloc( ref_face->f2n,
					  4 * ref_face_max(ref_face) * 
					  sizeof(REF_INT) );
      RNS(ref_face->f2n,"realloc ref_face->f2n NULL");
    }

  face = ref_face_n(ref_face);
  ref_face_n(ref_face)++;

  for(node=0;node<4;node++)
    ref_face_f2n(ref_face,node,face) = nodes[node];

  for(node=0;node<3;node++)
    RSS( ref_adj_add(ref_face_adj(ref_face),
		     ref_face_f2n(ref_face,node,face),
		     face ), "reg face node");

  if ( ref_face_f2n(ref_face,0,face) !=
       ref_face_f2n(ref_face,3,face) )
    RSS( ref_adj_add(ref_face_adj(ref_face),
		     ref_face_f2n(ref_face,3,face),
		     face ), "reg face node");

  return REF_SUCCESS;
}
