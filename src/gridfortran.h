
/* Michael A. Park
 * Computational Modeling & Simulation Branch
 * NASA Langley Research Center
 * Phone:(757)864-6604
 * Email:m.a.park@larc.nasa.gov
 */

#ifndef GRIDFORTRAN_H
#define GRIDFORTRAN_H

#include "refine_defs.h"

BEGIN_C_DECLORATION

void gridapiversion_( int *refine_api_version );

void gridcreate_( int *partId, int *nnode, double *x, double *y, double *z );
void gridfree_( void );

void gridinsertcells_( int *nodes_per_cell, int *ncell, int *c2n );
void gridinsertbc_( int *faceId, int *nodes_per_face, int *nface, int *f2n );
void gridsetmap_( int *nnode, double* map );
void gridsetimesh_( int *nnode, int* imesh );
void gridsetnodelocal2global_( int *partId, int *nnodeg, 
			      int *nnode, int *nnode0, int *local2global );
void gridsetnodepart_( int *nnode, int *part );
void gridfreezenode_( int *node );
void gridparallelloadcapri_( char *url, char *modeler, char *capriProject,
                             int *status );
void gridprojectallfaces_( void );
void gridtestcadparameters_( void );
void gridminar_( double *aspectratio );
void gridwritetecplotsurfacezone_( void );
void gridexportfast_( void );

void gridsetcostconstraint_( int *cost_constraint );
void gridconstrainsurfacenode_( void );

void gridparallelswap_( int *processor, double *ARlimit );
void gridparallelsmooth_( int *processor,
			  double *optimizationLimit, double *laplacianLimit,
                          int *geometryAllowed );
void gridparallelrelaxneg_( int *processor, int *geometryAllowed );
void gridparallelrelaxsurf_( int *processor );
void gridparalleladapt_( int *processor, 
			 double *minLength, double *maxLength );
void gridparallelpreproject_( int *processor );

void queuedumpsize_( int *nInt, int *nDouble );
void queuedump_( int *nInt, int *nDouble, int *ints, double *doubles );
void gridapplyqueue_( int *nInt, int *nDouble, int *ints, double *doubles );

void gridglobalnnode_( int *nnodeg );
void gridglobalshift_( int *oldnnodeg, int *newnnodeg, int *nodeoffset );
void gridrenumberglobalnodes_( int *nnode, int *new2old );

void gridnunusednodeglobal_( int *nunused );
void gridgetunusednodeglobal_( int *nunused, int *unused );
void gridjoinunusednodeglobal_( int *nunused, int *unused );
void gridcopyunusednodeglobal_( int *nunused, int *unused );
void grideliminateunusednodeglobal_( void );

void gridsortfun3d_( int *nnodes0, int *nnodes01, int *nnodesg );

void gridgetnodes_( int *nnode, int *l2g, double *x, double *y, double *z);
void gridgetimesh_( int *nnode, int *imesh);
void gridgetncell_( int *nodes_per_cell, int *ncell );
void gridgetcell_( int *nodes_per_cell, int *cell, int *nodes );
void gridgetbcsize_( int *ibound, int *nodes_per_face, int *nface );
void gridgetbc_( int *ibound, int *nodes_per_face, int *face, int *f2n );

void gridsetnaux_( int *naux );
void gridsetauxvector_( int *nnode, int *offset, double *x );
void gridsetauxmatrix_( int *ndim, int *nnode, int *offset, double *x );
void gridsetauxmatrix3_( int *ndim, int *nnode, int *offset, double *x );
void gridgetauxvector_( int *nnode, int *offset, double *x );
void gridgetauxmatrix_( int *ndim, int *nnode, int *offset, double *x );
void gridgetauxmatrix3_( int *ndim, int *nnode, int *offset, double *x );

void gridghostcount_( int *nproc, int *count );

void gridloadghostnodes_( int *nproc, int *clientindex,
			 int *clientsize, int *localnode, int *globalnode );
void gridloadglobalnodedata_( int *ndim, int *nnode, int *nodes, double *data );
void gridloadlocalnodes_( int *nnode, int *global, int *local );
void gridsetlocalnodedata_( int *ndim, int *nnode, int *nodes, double *data );

void gridmovesetprojectiondisp_( void );
void gridmoverelaxstartup_( int *relaxationScheme );
void gridmoverelaxstartstep_( double *position); 
void gridmoverelaxsubiter_( double *residual);
void gridmoverelaxshutdown_( void );
void gridmoveapplydisplacements_( void );

void gridmovedataleadingdim_( int *ndim );
void gridmoveinitializempitest_( void );
void gridmovecompletempitest_( void );
void gridmoveloadlocalnodedata_( int *ndim, int *nnode, 
				 int *nodes, double *data );
void gridmovesetlocalnodedata_( int *ndim, int *nnode, 
				int *nodes, double *data );

void gridmovefree_( void );

void gridgeomsize_( int *nGeomNode, int *nGeomEdge, int *nGeomFace );
void gridlocalboundnode_( int *nBoundNode );
void gridgeomedgeendpoints_( int *edgeId, int *endPoints );
void gridmaxedge_( int *maxedge );
void gridedge_( int *edge, int *edgeId, 
		int *globalnodes, int *nodeparts, 
		double *t, double *xyz);
void gridupdateedgegrid_(int *edgeId, int *nCurveNode, double *xyz, double *t);
void gridmaxface_( int *maxface );
void gridface_( int *face, int *faceId, 
		int *globalnodes, int *nodeparts, 
		double *uv, double *xyz);
void gridfaceedgecount_( int *faceId, int *faceEdgeCount );
void gridfaceedgel2g_( int *faceId, int *faceEdgeCount, int *local2global );

void gridupdategeometryface_( int *faceId, int *nnode, double *xyz, double *uv,
			      int *nface, int *f2n );

void gridcreateshellfromfaces_( void );

END_C_DECLORATION

#endif /* GRIDFORTRAN_H */
