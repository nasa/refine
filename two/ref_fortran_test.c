#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_fortran.h"

#include "ref_grid.h"
#include "ref_export.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_cell.h"
#include "ref_adj.h"

#include "ref_sort.h"
#include "ref_dict.h"
#include "ref_mpi.h"

int main( void )
{
  REF_INT nnodes;
  REF_INT nnodesg;
  REF_INT *l2g;
  REF_INT *part; 
  REF_INT partition;
  REF_DBL *x;
  REF_DBL *y; 
  REF_DBL *z;
  REF_INT *c2n;
  REF_INT *f2n;
  REF_DBL *m;

  REF_INT node_per_cell;
  REF_INT ncell;
  REF_INT node_per_face;
  REF_INT nface;
  REF_INT ibound;

  REF_INT node;

  nnodes = 4;
  nnodesg = 4;

  l2g  = (REF_INT *) malloc( sizeof(REF_INT) * nnodes );
  part = (REF_INT *) malloc( sizeof(REF_INT) * nnodes );
  x = (REF_DBL *) malloc( sizeof(REF_DBL) * nnodes );
  y = (REF_DBL *) malloc( sizeof(REF_DBL) * nnodes );
  z = (REF_DBL *) malloc( sizeof(REF_DBL) * nnodes );
  m = (REF_DBL *) malloc( sizeof(REF_DBL) * 6 * nnodes );

  l2g[0] = 3; part[0] = 1; x[0] = 0; y[0] = 1; z[0] = 0;
  l2g[1] = 4; part[1] = 1; x[1] = 0; y[1] = 0; z[1] = 1;
  l2g[2] = 1; part[2] = 2; x[2] = 0; y[2] = 0; z[2] = 0;
  l2g[3] = 2; part[3] = 2; x[3] = 1; y[3] = 0; z[3] = 0;

  for (node = 0; node < nnodes; node++)
    {
      m[0+6*node] = 1.0;
      m[1+6*node] = 0.0;
      m[2+6*node] = 0.0;
      m[3+6*node] = 1.0;
      m[4+6*node] = 0.0;
      m[5+6*node] = 1.0;
    }

  partition = 0;

  RSS(ref_init_node__(&nnodes, &nnodesg,
		      l2g, part, &partition,
		      x, y, z),"init node");

  node_per_cell = 4;
  ncell = 1;
  c2n = (REF_INT *) malloc( sizeof(REF_INT) * node_per_cell * ncell );

  c2n[0] = 1;
  c2n[1] = 2;
  c2n[2] = 3;
  c2n[3] = 4;

  RSS(ref_import_cell__( &node_per_cell, &ncell, c2n ),"import cell");

  node_per_face = 3;
  nface = 1;
  f2n = (REF_INT *) malloc( sizeof(REF_INT) * node_per_face * nface );
  f2n[0] = 3;
  f2n[1] = 4;
  f2n[2] = 1;

  ibound=1;
  RSS(ref_import_boundary__( &node_per_face, &nface, f2n, &ibound ),
      "import face");

  RSS(ref_import_metric__(&nnodes, m),"import metric");

  RSS(ref_free__(),"free");

  free(f2n);
  free(c2n);
  free(m);
  free(z);
  free(y);
  free(x);
  free(part);
  free(l2g);
  return 0;
}
