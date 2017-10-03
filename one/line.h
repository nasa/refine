
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#ifndef LINE_H
#define LINE_H

#include "refine_defs.h"

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
