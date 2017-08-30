
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:mike.park@nasa.gov
 */

#ifndef UNRAVEL_FORTRAN_H
#define UNRAVEL_FORTRAN_H

#include "refine_defs.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

BEGIN_C_DECLORATION
void FC_FUNC_(unravel_start,UNRAVEL_START)( int *unravel_api_version, int *status );
void FC_FUNC_(unravel_tet,UNRAVEL_TET)( int *c2n, double *x, double *y, double *z, int *status );
void FC_FUNC_(unravel_thaw,UNRAVEL_THAW)( int *nodeid, int *status );
void FC_FUNC_(unravel_it,UNRAVEL_IT)( int *status );
void FC_FUNC_(unravel_xyz,UNRAVEL_XYZ)( int *nodeid, double *x, double *y, double *z, int *status );
void FC_FUNC_(unravel_cleanup,UNRAVEL_CLEANUP)( int *status );
END_C_DECLORATION

#endif /* UNRAVEL_FORTRAN_H */

