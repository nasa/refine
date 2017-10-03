
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

#ifndef SORT_H
#define SORT_H

#include "refine_defs.h"

BEGIN_C_DECLORATION

void sortHeap( int length, int *arrayInput, int *sortedIndex );
void sortDoubleHeap( int length, double *arrayInput, int *sortedIndex );
int sortSearch( int length, int *arrayInput, int index );

END_C_DECLORATION

#endif /* QUEUE_H */
