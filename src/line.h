
/* a line in a mesh stored as node indices
 * 
 * Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov 
 */
  
/* $Id$ */

#ifndef LINE_H
#define LINE_H

#include "master_header.h"

BEGIN_C_DECLORATION

typedef struct Line Line;

struct Line {
  int length;
  int maxlength;
  int *nodes;
};

Line *lineCreate( void );
Line *lineInit( Line * );
void lineFree( Line * );
int lineLength( Line * );
Line *lineAddNode( Line *, int node );
int lineNode( Line *, int index );

typedef struct Lines Lines;

struct Lines {
  int number;
  int max;
  Line *line;
};

Lines *linesCreate( );
void linesFree( Lines * );
int linesNumber( Lines * );
Lines *linesAddNode( Lines *, int line, int node );
int linesNode( Lines *, int line, int index );
Lines *linesRenumber( Lines *, int *o2n );

Lines *linesSave( Lines *, char *filename );
Lines *linesLoad( Lines *, char *filename );

END_C_DECLORATION

#endif /* LINE_H */
