
/* Interp, a list items ranked in priorty
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

#include <stdlib.h>
#ifndef __APPLE__       /* Not needed on Mac OS X */
#include <malloc.h>
#endif
#include "interp.h"
#include "sort.h"

Interp* interpCreate( int function_id, int order )
{
  Interp *interp;

  interp = (Interp *)malloc( sizeof(Interp) );

  interp->function_id = function_id;
  interp->order = order;

  return interp;
}

void interpFree( Interp *interp )
{
  free( interp );
}

GridBool interpFunction( Interp *interp, double *xyz, double *func )
{
  (*func) = xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2];
  return TRUE;
}
