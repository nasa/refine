
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

REF_STATUS ref_html_create( REF_HTML *ref_html_ptr, char *filename );
REF_STATUS ref_html_free( REF_HTML ref_html );

#define ref_html_file(ref_html) ((ref_html)->file)

END_C_DECLORATION

#endif /* REF_HTML_H */
