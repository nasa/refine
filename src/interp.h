
/* Interp, a list items ranked in priorty
 *
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

/* $Id: interp.h 2867 2005-07-18 17:32:34Z mikepark $ */

#ifndef INTERP_H
#define INTERP_H

#include "refine_defs.h"
#include <stdio.h>

BEGIN_C_DECLORATION

typedef struct Interp Interp;

struct Interp {
  int function_id;
};

Interp *interpCreate( int function_id );
void interpFree( Interp * );

GridBool interpFunction( Interp *, double *xyz, double *func );

END_C_DECLORATION

#endif /* INTERP_H */
