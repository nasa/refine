
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_EGADS
#include "egads.h"
#endif

#include "ref_geom.h"

#include "ref_dict.h"
#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_malloc.h"
#include "ref_adapt.h"
#include "ref_matrix.h"

REF_STATUS ref_geom_egads_fixture( char *filename )
{
#ifdef HAVE_EGADS
  ego context;
  int stype;
  double data[6];
  ego body;
  ego model = NULL;
  
  REIS( EGADS_SUCCESS, EG_open(&context), "EG open");
  stype = 1; /* box */
  data[0]=0.0; /* corner */
  data[1]=0.0;
  data[2]=0.0;
  data[3]=1.0; /* length of sides */
  data[4]=1.0;
  data[5]=1.0;
  REIS( EGADS_SUCCESS,
	EG_makeSolidBody(context, stype, data, &body), "EG box");
  REIS( EGADS_SUCCESS,
	EG_makeTopology(context, NULL, MODEL, 0,
			NULL, 1, &body,
			NULL, &model), "topo");
  REIS( EGADS_SUCCESS,
	EG_saveModel(model, filename), "save");
  printf("wrote EGADS project %s\n",filename);

#else
  printf("No EGADS linked for %s\n",filename);
#endif
  
  return REF_SUCCESS;
}

REF_STATUS ref_geom_from_egads( REF_GRID *ref_grid_ptr, char *filename )
{
#ifdef HAVE_EGADS
  REF_GRID ref_grid;
  REF_NODE ref_node;
  ego context;
  ego model = NULL;
  double params[3], box[6], size;
  ego geom, *bodies, *dum;
  int oclass, mtype, nbody, *senses, j;
  ego solid, tess;
  int tess_status, nvert;
  
  printf("EGAGS project %s\n",filename);
  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  REIS( EGADS_SUCCESS, EG_open(&context), "EG open");
  REIS( EGADS_SUCCESS, EG_loadModel(context, 0, filename, &model), "EG load");
  REIS( EGADS_SUCCESS, EG_getBoundingBox(model, box), "EG bounding box");
  size = sqrt((box[0]-box[3])*(box[0]-box[3]) +
	      (box[1]-box[4])*(box[1]-box[4]) +
	      (box[2]-box[5])*(box[2]-box[5]));
  printf(" box %f %f %f %f %f %f\n",box[0],box[1],box[2],box[3],box[4],box[5]);

  params[0] =  0.25*size; /*spacing*/
  params[1] =  0.001*size;
  params[2] = 15.0;

  REIS( EGADS_SUCCESS,
	EG_getTopology(model, &geom, &oclass, &mtype, NULL,
		       &nbody, &bodies, &senses), "EG topo bodies");
  printf(" %d bodies\n",nbody);
  REIS( EGADS_SUCCESS,
	EG_getTopology(bodies[0], &geom, &oclass, &mtype,
		       NULL, &j, &dum, &senses), "EG topo body type");
  REIS( SOLIDBODY, mtype, "expected SOLIDBODY" );
  REIS( EGADS_SUCCESS,
	EG_makeTopology(context, NULL, BODY, SOLIDBODY, NULL, j, dum,
			NULL, &solid), "EG topo solid");
  
  REIS( EGADS_SUCCESS,
	EG_makeTessBody(solid, params, &tess), "EG tess");

  REIS( EGADS_SUCCESS,
	EG_statusTessBody(tess, &geom, &tess_status, &nvert), "EG tess");
  REIS( 1, tess_status, "tess not closed" );
  printf(" %d global vertex\n",nvert);
  
  
  
#else
  printf("returning empty grid, No EGADS linked for %s\n",filename);
  RSS( ref_grid_create( ref_grid_ptr ), "create grid");  
#endif
  
  return REF_SUCCESS;
}
