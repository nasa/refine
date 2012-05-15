
#include <stdlib.h>
#include <stdio.h>

#include <string.h>

#include "ref_plot3d.h"
#include "ref_malloc.h"

REF_STATUS ref_plot3d_from_file( REF_PLOT3D *ref_plot3d_ptr, char *filename )
{
  REF_PLOT3D ref_plot3d;
  FILE *file;
  REF_INT fortran_record_size;

  ref_malloc( *ref_plot3d_ptr, 1, REF_PLOT3D_STRUCT );
  ref_plot3d = *ref_plot3d_ptr;

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "nnode" );
  REIS( 1*4, fortran_record_size, "header start record size" );

  RES( 1, fread( &ref_plot3d_ngrid(ref_plot3d), sizeof(REF_INT), 1,file),"ng" );

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "nnode" );
  REIS( 1*4, fortran_record_size, "header end record size" );

  fclose(file);

  return REF_SUCCESS;
}
