
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

int gridcreate_( int *nnode, double *x, double *y, double *z,
		 int *ncell, int *maxcell, int *c2n );
int gridfree_( );
int gridinsertboundary_( int *faceId, int *nnode, int *nodedim, int *inode, 
			 int *nface, int *dim1, int *dim2, int *f2n );
int gridsetmap_( int *nnode, double* map );
int gridsetnodelocal2global_( int *partId, int *nnodeg, 
			      int *nnode, int *nnode0, int *local2global );
int gridsetnodepart_( int *nnode, int *part );
int gridsetcelllocal2global_( int *ncell, int *local2global );
int gridswap_( );
int gridsmoothvolume_( );
int gridadaptwithoutcad_( double *minLength, double *maxLength );
int gridwritetecplotsurfacezone_( );

int gridparalleladaptwithoutcad_( int *processor, 
				  double *minLength, double *maxLength );
int queuedumpsize_( int *nInt, int *nDouble );
int queuedump_( int *nInt, int *nDouble, int *ints, double *doubles );
int gridapplyqueue_( int *nInt, int *nDouble, int *ints, double *doubles );

int gridsize_( int *nnodeg, int *ncellg );

int gridglobalshift_( int *oldnnodeg, int *newnnodeg, int *nodeoffset,
		      int *oldncellg, int *newncellg, int *celloffset );

int gridnunusedcellglobal_( int *nunused );
int gridgetunusedcellglobal_( int *nunused, int *unused );
int gridjoinunusedcellglobal_( int *nunused, int *unused );

END_C_DECLORATION

#endif /* GRIDFORTRAN_H */
