
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id$ */

#ifndef FORTRAN_H
#define FORTRAN_H

#include "master_header.h"

BEGIN_C_DECLORATION

END_C_DECLORATION

#endif /* FORTRAN_H */

int gridcreate( int *maxnode, int *maxcell, int *maxface, int *maxedge );
int gridcreate_( int *maxnode, int *maxcell, int *maxface, int *maxedge );
int gridcreate__( int *maxnode, int *maxcell, int *maxface, int *maxedge );
int GRIDCREATE( int *maxnode, int *maxcell, int *maxface, int *maxedge );
