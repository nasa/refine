
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
void unravel_start_( int *unravel_api_version, int *status );
void unravel_tet_( int *c2n, double *x, double *y, double *z, int *status );
void unravel_thaw_( int *nodeid, int *status );
void unravel_it_( int *status );
void unravel_xyz_( int *nodeid, double *x, double *y, double *z, int *status );
void unravel_cleanup_( int *status );
END_C_DECLORATION

#endif /* UNRAVEL_FORTRAN_H */

