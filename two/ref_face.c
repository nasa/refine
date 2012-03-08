
#include <stdlib.h>
#include <stdio.h>

#include "ref_face.h"
#include "ref_sort.h"
#include "ref_math.h"

#include "ref_malloc.h"

REF_STATUS ref_face_create( REF_FACE *ref_face_ptr, REF_GRID ref_grid )
{
  REF_FACE ref_face;
  REF_INT group, cell, cell_face;
  REF_INT node;
  REF_INT nodes[4];
  REF_CELL ref_cell;

  ref_malloc( *ref_face_ptr, 1, REF_FACE_STRUCT );

  ref_face = *ref_face_ptr;

  ref_face_n(ref_face) = 0;
  ref_face_max(ref_face) = 10;

  ref_malloc( ref_face->f2n, 4 * ref_face_max(ref_face), REF_INT );

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
  ref_free( ref_face->f2n );

  ref_free( ref_face );

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
      RSS( ref_sort_insertion( 4, orig, sort ), "sort" );
      printf("sort %4d : %4d %4d %4d %4d\n",face,
	     sort[0],sort[1],sort[2],sort[3]);
    }
  return REF_SUCCESS;
}

static REF_STATUS ref_face_make_canonical( REF_INT *original, 
					   REF_INT *canonical )
{
  RSS( ref_sort_insertion( 4, original, canonical ), "sort" );

  if ( canonical[0] == canonical[1] )
    {
      canonical[0] = REF_EMPTY;
      return REF_SUCCESS;
    }

  if ( canonical[1] == canonical[2] )
    {
      canonical[1] = canonical[0];
      canonical[0] = REF_EMPTY;
      return REF_SUCCESS;
    }

  if ( canonical[2] == canonical[3] )
    {
      canonical[2] = canonical[1];
      canonical[1] = canonical[0];
      canonical[0] = REF_EMPTY;
      return REF_SUCCESS;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_face_with( REF_FACE ref_face, REF_INT *nodes, REF_INT *face )
{
  REF_INT item, ref, node;
  REF_INT target[4], canidate[4], orig[4];

  (*face) = REF_EMPTY;

  RSS( ref_face_make_canonical( nodes, target ), "canonical" );

  each_ref_adj_node_item_with_ref( ref_face_adj(ref_face), nodes[0], item, ref)
    {
      for(node=0;node<4;node++)
	orig[node]=ref_face_f2n(ref_face,node,ref);
      RSS( ref_face_make_canonical( orig, canidate ), "canonical" );
      if ( target[0] == canidate[0] &&
	   target[1] == canidate[1] &&
	   target[2] == canidate[2] &&
	   target[3] == canidate[3] )
	{
	  (*face) = ref;
	  return REF_SUCCESS;
	}
    }

  return REF_NOT_FOUND;
}

REF_STATUS ref_face_spanning( REF_FACE ref_face, REF_INT node0, REF_INT node1, 
			      REF_INT *face )
{
  REF_INT item, ref;

  (*face) = REF_EMPTY;

  each_ref_adj_node_item_with_ref( ref_face_adj(ref_face), node0, item, ref)
    {
      if ( ( node0 == ref_face_f2n(ref_face,0,ref) ||
	     node0 == ref_face_f2n(ref_face,1,ref) ||
	     node0 == ref_face_f2n(ref_face,2,ref) ||
	     node0 == ref_face_f2n(ref_face,3,ref) ) &&
	   ( node1 == ref_face_f2n(ref_face,0,ref) ||
	     node1 == ref_face_f2n(ref_face,1,ref) ||
	     node1 == ref_face_f2n(ref_face,2,ref) ||
	     node1 == ref_face_f2n(ref_face,3,ref) ) )
	{
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
      ref_realloc( ref_face->f2n, 4 * ref_face_max(ref_face), REF_INT );
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

REF_STATUS ref_face_normal( REF_DBL *xyz0, REF_DBL *xyz1, 
			    REF_DBL *xyz2, REF_DBL *xyz3, 
			    REF_DBL *normal )
{
  REF_DBL edge1[3], edge2[3];

  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 0.0;

  edge1[0] = xyz2[0] - xyz1[0];
  edge1[1] = xyz2[1] - xyz1[1];
  edge1[2] = xyz2[2] - xyz1[2];

  edge2[0] = xyz0[0] - xyz1[0];
  edge2[1] = xyz0[1] - xyz1[1];
  edge2[2] = xyz0[2] - xyz1[2];

  normal[0] += 0.5*(edge1[1]*edge2[2] - edge1[2]*edge2[1]);
  normal[1] += 0.5*(edge1[2]*edge2[0] - edge1[0]*edge2[2]);
  normal[2] += 0.5*(edge1[0]*edge2[1] - edge1[1]*edge2[0]);

  edge1[0] = xyz0[0] - xyz3[0];
  edge1[1] = xyz0[1] - xyz3[1];
  edge1[2] = xyz0[2] - xyz3[2];

  edge2[0] = xyz2[0] - xyz3[0];
  edge2[1] = xyz2[1] - xyz3[1];
  edge2[2] = xyz2[2] - xyz3[2];

  normal[0] += 0.5*(edge1[1]*edge2[2] - edge1[2]*edge2[1]);
  normal[1] += 0.5*(edge1[2]*edge2[0] - edge1[0]*edge2[2]);
  normal[2] += 0.5*(edge1[0]*edge2[1] - edge1[1]*edge2[0]);

  return REF_SUCCESS;

}

REF_STATUS ref_face_open_node( REF_DBL *xyz0, REF_DBL *xyz1, 
			       REF_DBL *xyz2, REF_DBL *xyz3, 
			       REF_INT *open_node )
{
  REF_DBL edge1[3], edge2[3];
  REF_DBL open_dot;

  edge1[0] = xyz0[0] - xyz3[0];
  edge1[1] = xyz0[1] - xyz3[1];
  edge1[2] = xyz0[2] - xyz3[2];

  edge2[0] = xyz1[0] - xyz0[0];
  edge2[1] = xyz1[1] - xyz0[1];
  edge2[2] = xyz1[2] - xyz0[2];

  open_dot = ref_math_dot( edge2, edge1 );
  *open_node = 0;

  edge1[0] = xyz1[0] - xyz0[0];
  edge1[1] = xyz1[1] - xyz0[1];
  edge1[2] = xyz1[2] - xyz0[2];

  edge2[0] = xyz2[0] - xyz1[0];
  edge2[1] = xyz2[1] - xyz1[1];
  edge2[2] = xyz2[2] - xyz1[2];

  if ( ref_math_dot( edge2, edge1 ) > open_dot )
    {
      open_dot = ref_math_dot( edge2, edge1 );
      *open_node = 1;
    }

  edge1[0] = xyz2[0] - xyz1[0];
  edge1[1] = xyz2[1] - xyz1[1];
  edge1[2] = xyz2[2] - xyz1[2];

  edge2[0] = xyz3[0] - xyz2[0];
  edge2[1] = xyz3[1] - xyz2[1];
  edge2[2] = xyz3[2] - xyz2[2];

  if ( ref_math_dot( edge2, edge1 ) > open_dot )
    {
      open_dot = ref_math_dot( edge2, edge1 );
      *open_node = 2;
    }

  edge1[0] = xyz3[0] - xyz2[0];
  edge1[1] = xyz3[1] - xyz2[1];
  edge1[2] = xyz3[2] - xyz2[2];

  edge2[0] = xyz0[0] - xyz3[0];
  edge2[1] = xyz0[1] - xyz3[1];
  edge2[2] = xyz0[2] - xyz3[2];

  if ( ref_math_dot( edge2, edge1 ) > open_dot )
    {
      open_dot = ref_math_dot( edge2, edge1 );
      *open_node = 3;
    }

  return REF_SUCCESS;

}
