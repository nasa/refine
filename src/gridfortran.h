
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

void gridcreate_( int *partId, int *nnode, double *x, double *y, double *z,
		 int *ncell, int *maxcell, int *c2n );
void gridfree_( );
void gridinsertboundary_( int *faceId, int *nnode, int *nodedim, int *inode, 
			 int *nface, int *dim1, int *dim2, int *f2n );
void gridsetmap_( int *nnode, double* map );
void gridsetnodelocal2global_( int *partId, int *nnodeg, 
			      int *nnode, int *nnode0, int *local2global );
void gridsetnodepart_( int *nnode, int *part );
void gridsetcelllocal2global_( int *ncellg, int *ncell, int *local2global );
void gridprojectallfaces_( );
void gridswap_( );
void gridsmoothvolume_( );
void gridsmoothfaceinterior_( int *processor );
void gridadaptwithoutcad_( double *minLength, double *maxLength );
void gridwritetecplotsurfacezone_( );
void gridexportfast_( );

void gridparalleladaptwithoutcad_( int *processor, 
				  double *minLength, double *maxLength );
void gridparallelswap_( int *processor );
void queuedumpsize_( int *nInt, int *nDouble );
void queuedump_( int *nInt, int *nDouble, int *ints, double *doubles );
void gridapplyqueue_( int *nInt, int *nDouble, int *ints, double *doubles );

void gridsize_( int *nnodeg, int *ncellg );

void gridglobalshift_( int *oldnnodeg, int *newnnodeg, int *nodeoffset,
		      int *oldncellg, int *newncellg, int *celloffset );

void gridnunusedcellglobal_( int *nunused );
void gridgetunusedcellglobal_( int *nunused, int *unused );
void gridjoinunusedcellglobal_( int *nunused, int *unused );
void grideliminateunusedcellglobal_( );

void gridsortfun3d_( int *nnodes0, int *nnodes01, int *nnodesg, 
		    int *ncell, int *ncellg );
void gridgetnodes_( int *nnode, int *l2g, double *x, double *y, double *z);
void gridgetcell_( int *cell, int *nodes, int *global );
void gridgetbcsize_( int *ibound, int *nface );
void gridgetbc_( int *ibound, int *nface, int *ndim, int *f2n );

void gridsetnaux_( int *naux );
void gridsetauxvector_( int *nnode, int *offset, double *x );
void gridsetauxmatrix_( int *ndim, int *nnode, int *offset, double *x );
void gridsetauxmatrix3_( int *ndim, int *nnode, int *offset, double *x );
void gridgetauxvector_( int *nnode, int *offset, double *x );
void gridgetauxmatrix_( int *ndim, int *nnode, int *offset, double *x );
void gridgetauxmatrix3_( int *ndim, int *nnode, int *offset, double *x );

void gridghostcount_( int *nproc, int *count );

void gridloadghostnodes_( int *nproc, int *clientindex,
			 int *clientsize, int *localnode, int *globalnode );
void gridloadglobalnodedata_( int *ndim, int *nnode, int *nodes, double *data );
void gridsetlocalnodedata_( int *ndim, int *nnode, int *nodes, double *data );

END_C_DECLORATION

#endif /* GRIDFORTRAN_H */
