
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#ifndef GRIDFORTRAN_H
#define GRIDFORTRAN_H

#include "master_header.h"

BEGIN_C_DECLORATION

int gridcreate_( int *partId, int *nnode, double *x, double *y, double *z,
		 int *ncell, int *maxcell, int *c2n );
int gridfree_( );
int gridinsertboundary_( int *faceId, int *nnode, int *nodedim, int *inode, 
			 int *nface, int *dim1, int *dim2, int *f2n );
int gridsetmap_( int *nnode, double* map );
int gridsetnodelocal2global_( int *partId, int *nnodeg, 
			      int *nnode, int *nnode0, int *local2global );
int gridsetnodepart_( int *nnode, int *part );
int gridsetcelllocal2global_( int *ncellg, int *ncell, int *local2global );
int gridswap_( );
int gridsmoothvolume_( );
int gridadaptwithoutcad_( double *minLength, double *maxLength );
int gridwritetecplotsurfacezone_( );

int gridparalleladaptwithoutcad_( int *processor, 
				  double *minLength, double *maxLength );
int gridparallelswap_( int *processor );
int queuedumpsize_( int *nInt, int *nDouble );
int queuedump_( int *nInt, int *nDouble, int *ints, double *doubles );
int gridapplyqueue_( int *nInt, int *nDouble, int *ints, double *doubles );

int gridsize_( int *nnodeg, int *ncellg );

int gridglobalshift_( int *oldnnodeg, int *newnnodeg, int *nodeoffset,
		      int *oldncellg, int *newncellg, int *celloffset );

int gridnunusedcellglobal_( int *nunused );
int gridgetunusedcellglobal_( int *nunused, int *unused );
int gridjoinunusedcellglobal_( int *nunused, int *unused );
int grideliminateunusedcellglobal_( );

int gridsortfun3d_( int *nnodes0, int *nnodes01, int *nnodesg, 
		    int *ncell, int *ncellg );
int gridgetnodes_( int *nnode, int *l2g, double *x, double *y, double *z);
int gridgetcell_( int *cell, int *nodes, int *global );
int gridgetbcsize_( int *ibound, int *nface );
int gridgetbc_( int *ibound, int *nface, int *ndim, int *f2n );

int gridsetnaux_( int *naux );
int gridsetauxvector_( int *nnode, int *offset, double *x );
int gridsetauxmatrix_( int *ndim, int *nnode, int *offset, double *x );
int gridsetauxmatrix3_( int *ndim, int *nnode, int *offset, double *x );
int gridgetauxvector_( int *nnode, int *offset, double *x );
int gridgetauxmatrix_( int *ndim, int *nnode, int *offset, double *x );
int gridgetauxmatrix3_( int *ndim, int *nnode, int *offset, double *x );

int gridghostcount_( int *nproc, int *count );

END_C_DECLORATION

#endif /* GRIDFORTRAN_H */
