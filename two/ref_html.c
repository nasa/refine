
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

  fprintf(f,"<html>\n");
  fprintf(f,"  <head>\n");
  fprintf(f,"    <title>refine export</title>\n");
  fprintf(f,"    <link rel='stylesheet' type='text/css'\n");
  fprintf(f,"          href='http://www.x3dom.org/download/x3dom.css'>\n");
  fprintf(f,"    </link>\n");
  fprintf(f,"    <script type='text/javascript'\n");
  fprintf(f,"          src='http://www.x3dom.org/download/x3dom.js'>\n");
  fprintf(f,"    </script>\n");
  fprintf(f,"    <style>\n");
  fprintf(f,"      x3d {width:100%%;height:100%%;border:none}\n");
  fprintf(f,"      body {margin:0;width:100%%;height:100%%;}\n");
  fprintf(f,"    </style>\n");
  fprintf(f,"  </head>\n");
  fprintf(f,"  <body id='body'>\n");
  fprintf(f,"    <a href=\"http://x3dom.org/docs/dev/navigation.html\">\n");
  fprintf(f,"camera control help\n");
  fprintf(f,"    </a>\n");
  fprintf(f,"    <x3d id='x3d'><scene><shape>\n");

  return REF_SUCCESS;
}

REF_STATUS ref_html_free( REF_HTML ref_html )
{
  if ( NULL != (void *)ref_html_file(ref_html) ) 
    {
      fprintf(ref_html_file(ref_html),"    </shape></scene></x3d>\n");
      fprintf(ref_html_file(ref_html),"  </body>\n");
      fprintf(ref_html_file(ref_html),"</html>\n");
      fclose( ref_html_file(ref_html) );
      ref_html_file(ref_html) = NULL;
    }
  ref_free( ref_html );
  return REF_SUCCESS;
}
