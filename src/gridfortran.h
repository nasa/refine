
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
int gridswap_( );
int gridsmoothvolume_( );
int gridadaptwithoutcad_( double *minLength, double *maxLength );
int gridwritetecplotsurfacezone_( );

END_C_DECLORATION

#endif /* GRIDFORTRAN_H */
