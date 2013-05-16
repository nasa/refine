
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_html.h"
#include "ref_malloc.h"

REF_STATUS ref_html_create( REF_HTML *ref_html_ptr, char *filename )
{
  FILE *f;

  f = fopen(filename,"w");
  if (NULL == (void *)f) printf("unable to open %s\n",filename);
  RNS(f, "unable to open file" );

  ref_malloc( *ref_html_ptr, 1, REF_HTML_STRUCT );
  ref_html_file(*ref_html_ptr) = f;

  return REF_SUCCESS;
}

REF_STATUS ref_html_free( REF_HTML ref_html )
{
  if ( NULL != (void *)ref_html_file(ref_html) ) 
    fclose( ref_html_file(ref_html) );
  ref_html_file(ref_html) = NULL;
  ref_free( ref_html );
  return REF_SUCCESS;
}
