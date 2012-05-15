
#include <stdlib.h>
#include <stdio.h>

#include <string.h>

#include "ref_plot3d.h"
#include "ref_malloc.h"

#include "ref_endian.h"

/* http://cfl3d.larc.nasa.gov/Cfl3dv6/V5Manual/FileForm.pdf */

REF_STATUS ref_plot3d_from_file( REF_PLOT3D *ref_plot3d_ptr, char *filename )
{
  REF_PLOT3D ref_plot3d;
  FILE *file;
  REF_INT fortran_record_size;
  REF_INT n;

  ref_malloc( *ref_plot3d_ptr, 1, REF_PLOT3D_STRUCT );
  ref_plot3d = *ref_plot3d_ptr;

  file = fopen(filename,"r");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "srt rec" );
  SWAP_INT(fortran_record_size);
  REIS( 1*4, fortran_record_size, "header start record size" );

  RES( 1, fread( &ref_plot3d_ngrid(ref_plot3d), sizeof(REF_INT), 1,file),"ng" );
  SWAP_INT(ref_plot3d_ngrid(ref_plot3d));

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "end rec" );
  SWAP_INT(fortran_record_size);
  REIS( 1*4, fortran_record_size, "header end record size" );

  ref_malloc( ref_plot3d->idim, ref_plot3d_ngrid(ref_plot3d), REF_INT );
  ref_malloc( ref_plot3d->jdim, ref_plot3d_ngrid(ref_plot3d), REF_INT );
  ref_malloc( ref_plot3d->kdim, ref_plot3d_ngrid(ref_plot3d), REF_INT );

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "srt rec" );
  SWAP_INT(fortran_record_size);
  REIS( ref_plot3d_ngrid(ref_plot3d)*3*4, fortran_record_size, 
	"dims start record size" );

  for (n=0;n<ref_plot3d_ngrid(ref_plot3d);n++)
    {
      RES( 1, fread( &(ref_plot3d->idim[n]), sizeof(REF_INT), 1,file),"id" );
      SWAP_INT(ref_plot3d->idim[n]);
      RES( 1, fread( &(ref_plot3d->jdim[n]), sizeof(REF_INT), 1,file),"jd" );
      SWAP_INT(ref_plot3d->jdim[n]);
      RES( 1, fread( &(ref_plot3d->kdim[n]), sizeof(REF_INT), 1,file),"kd" );
      SWAP_INT(ref_plot3d->kdim[n]);
    }

  RES( 1, fread( &fortran_record_size, sizeof(REF_INT), 1, file ), "end rec" );
  SWAP_INT(fortran_record_size);
  REIS( ref_plot3d_ngrid(ref_plot3d)*3*4, fortran_record_size, 
	"dims end record size" );

  fclose(file);

  return REF_SUCCESS;
}
