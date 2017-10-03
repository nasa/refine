
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

#ifndef REF_HTML_H
#define REF_HTML_H

#include <stdio.h>

#include "ref_defs.h"

BEGIN_C_DECLORATION
typedef struct REF_HTML_STRUCT REF_HTML_STRUCT;
typedef REF_HTML_STRUCT * REF_HTML;
END_C_DECLORATION

BEGIN_C_DECLORATION

struct REF_HTML_STRUCT {
  FILE *file;
};

REF_STATUS ref_html_create( REF_HTML *ref_html_ptr, const char *filename );
REF_STATUS ref_html_diagonal_system( REF_HTML ref_html, 
				     REF_DBL *origin,
				     REF_DBL *diagonal_system );
REF_STATUS ref_html_free( REF_HTML ref_html );

#define ref_html_file(ref_html) ((ref_html)->file)

END_C_DECLORATION

#endif /* REF_HTML_H */
