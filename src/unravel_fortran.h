
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:mike.park@nasa.gov
 */

#ifndef UNRAVEL_FORTRAN_H
#define UNRAVEL_FORTRAN_H

#include "refine_defs.h"

BEGIN_C_DECLORATION
void unravel_start_( void );
void unravel_tet_( int *c2n, double *x, double *y, double *z );
void unravel_thaw_( int *nodeid );
void unravel_it_( int *success );
void unravel_xyz_( int *nodeid, double *x, double *y, double *z );
void unravel_cleanup_( void );
END_C_DECLORATION

#endif /* UNRAVEL_FORTRAN_H */

