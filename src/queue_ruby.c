
#include "ruby.h"
#include "queue.h"

#define GET_QUEUE_FROM_SELF Queue *queue; Data_Get_Struct( self, Queue, queue );

static void queue_free( void *queue )
{
  queueFree( queue );
}

VALUE queue_new( VALUE class, VALUE nodeSize )
{
  Queue *queue;
  VALUE obj;
  queue = queueCreate( NUM2INT(nodeSize) );
  obj = Data_Wrap_Struct( class, 0, queue_free, queue );
  return obj;
}

VALUE queue_reset(VALUE self )
{
  GET_QUEUE_FROM_SELF;
  return (queue==queueReset(queue)?self:Qnil);
}

VALUE queue_nodeSize( VALUE self )
{
  GET_QUEUE_FROM_SELF;
  return INT2NUM(queueNodeSize(queue));
}

VALUE queue_transactions( VALUE self )
{
  GET_QUEUE_FROM_SELF;
  return INT2NUM(queueTransactions(queue));
}

VALUE queue_newTransaction( VALUE self )
{
  GET_QUEUE_FROM_SELF;
  return (queue==queueNewTransaction(queue)?self:Qnil);
}

/* ****************************** cells ****************************** */

VALUE queue_removeCell( VALUE self, VALUE rb_nodes, VALUE rb_nodeParts )
{
  int i, nodes[4], nodeParts[4];
  GET_QUEUE_FROM_SELF;
  for (i=0;i<4;i++) nodes[i]=NUM2INT(rb_ary_entry(rb_nodes,i));
  for (i=0;i<4;i++) nodeParts[i]=NUM2INT(rb_ary_entry(rb_nodeParts,i));
  return (queue==queueRemoveCell(queue,nodes,nodeParts)?self:Qnil);
}

VALUE queue_removedCells( VALUE self, VALUE transaction )
{
  GET_QUEUE_FROM_SELF;
  return INT2NUM(queueRemovedCells(queue,NUM2INT(transaction)));
}

VALUE queue_removedCellNodes( VALUE self, VALUE index )
{
  int i, nodes[4];
  VALUE rb_nodes;
  GET_QUEUE_FROM_SELF;
  if (queue != queueRemovedCellNodes(queue,NUM2INT(index),nodes)) return Qnil;
  rb_nodes = rb_ary_new2(4);
  for (i=0;i<4;i++) rb_ary_store(rb_nodes,i,INT2NUM(nodes[i]));
  return rb_nodes;
}

VALUE queue_removedCellNodeParts( VALUE self, VALUE index )
{
  int i, nodeParts[4];
  VALUE rb_nodeParts;
  GET_QUEUE_FROM_SELF;
  if (queue != queueRemovedCellNodeParts(queue,NUM2INT(index),nodeParts)) 
    return Qnil;
  rb_nodeParts = rb_ary_new2(4);
  for (i=0;i<4;i++) rb_ary_store(rb_nodeParts,i,INT2NUM(nodeParts[i]));
  return rb_nodeParts;
}

VALUE queue_addCell( VALUE self, VALUE rb_nodes, VALUE rb_cellId, 
		     VALUE rb_nodeParts, VALUE rb_xyzs )
{
  int i, nodes[4], cellId, nodeParts[4];
  double *xyzs;
  VALUE result;
  GET_QUEUE_FROM_SELF;
  for (i=0;i< 4;i++) {
    nodes[i]=NUM2INT(rb_ary_entry(rb_nodes,i));
    nodeParts[i]=NUM2INT(rb_ary_entry(rb_nodeParts,i));
  }
  cellId = NUM2INT(rb_cellId);
  xyzs = malloc(4*queueNodeSize(queue)*sizeof(double));
  for (i=0;i<4*queueNodeSize(queue);i++) 
    xyzs[i] = NUM2DBL(rb_ary_entry(rb_xyzs,i));
  result = (queue==queueAddCell(queue,nodes,cellId,nodeParts,xyzs)?self:Qnil);
  free(xyzs);
  return result;
}

VALUE queue_addedCells( VALUE self, VALUE transaction )
{
  GET_QUEUE_FROM_SELF;
  return INT2NUM(queueAddedCells(queue,NUM2INT(transaction)));
}

VALUE queue_addedCellNodes( VALUE self, VALUE index )
{
  int i, nodes[4];
  VALUE rb_nodes;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedCellNodes(queue,NUM2INT(index),nodes)) return Qnil;
  rb_nodes = rb_ary_new2(4);
  for (i=0;i<4;i++) rb_ary_store(rb_nodes,i,INT2NUM(nodes[i]));
  return rb_nodes;
}

VALUE queue_addedCellId( VALUE self, VALUE index )
{
  int cellId;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedCellId(queue,NUM2INT(index),&cellId)) return Qnil;
  return INT2NUM(cellId);
}

VALUE queue_addedCellNodeParts( VALUE self, VALUE index )
{
  int i, nodeParts[4];
  VALUE rb_nodeParts;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedCellNodeParts(queue,NUM2INT(index),nodeParts)) return Qnil;
  rb_nodeParts = rb_ary_new2(4);
  for (i=0;i<4;i++) rb_ary_store(rb_nodeParts,i,INT2NUM(nodeParts[i]));
  return rb_nodeParts;
}

VALUE queue_addedCellXYZs( VALUE self, VALUE index )
{
  int i;
  double *xyzs;
  VALUE rb_xyzs;
  GET_QUEUE_FROM_SELF;
  xyzs = malloc(4*queueNodeSize(queue)*sizeof(double));
  if (queue != queueAddedCellXYZs(queue,NUM2INT(index),xyzs)) {
    free (xyzs);
    return Qnil;
  }else{
    rb_xyzs = rb_ary_new2(4*queueNodeSize(queue));
    for (i=0;i<4*queueNodeSize(queue);i++) 
      rb_ary_store(rb_xyzs,i,rb_float_new(xyzs[i]));
    return rb_xyzs;
  }
}

VALUE queue_totalRemovedCells( VALUE self )
{
  GET_QUEUE_FROM_SELF;
  return INT2NUM(queueTotalRemovedCells(queue));
}

/* ****************************** faces ****************************** */

VALUE queue_removeFace( VALUE self, VALUE rb_nodes, VALUE rb_nodeParts )
{
  int i, nodes[3], nodeParts[3];
  GET_QUEUE_FROM_SELF;
  for (i=0;i<3;i++) nodes[i]=NUM2INT(rb_ary_entry(rb_nodes,i));
  for (i=0;i<3;i++) nodeParts[i]=NUM2INT(rb_ary_entry(rb_nodeParts,i));
  return (queue==queueRemoveFace(queue,nodes,nodeParts)?self:Qnil);
}

VALUE queue_removedFaces( VALUE self, VALUE transaction )
{
  GET_QUEUE_FROM_SELF;
  return INT2NUM(queueRemovedFaces(queue,NUM2INT(transaction)));
}

VALUE queue_removedFaceNodes( VALUE self, VALUE index )
{
  int i, nodes[3];
  VALUE rb_nodes;
  GET_QUEUE_FROM_SELF;
  if (queue != queueRemovedFaceNodes(queue,NUM2INT(index),nodes)) return Qnil;
  rb_nodes = rb_ary_new2(3);
  for (i=0;i<3;i++) rb_ary_store(rb_nodes,i,INT2NUM(nodes[i]));
  return rb_nodes;
}

VALUE queue_removedFaceNodeParts( VALUE self, VALUE index )
{
  int i, nodeParts[3];
  VALUE rb_nodeParts;
  GET_QUEUE_FROM_SELF;
  if (queue != queueRemovedFaceNodeParts(queue,NUM2INT(index),nodeParts)) 
    return Qnil;
  rb_nodeParts = rb_ary_new2(3);
  for (i=0;i<3;i++) rb_ary_store(rb_nodeParts,i,INT2NUM(nodeParts[i]));
  return rb_nodeParts;
}

VALUE queue_addFace( VALUE self, VALUE rb_nodes, VALUE rb_faceId, 
		     VALUE rb_nodeParts, VALUE rb_uvs )
{
  int i, nodes[3], faceId, nodeParts[3];
  double uvs[6];
  GET_QUEUE_FROM_SELF;
  for (i=0;i< 3;i++) {
    nodes[i]=NUM2INT(rb_ary_entry(rb_nodes,i));
    nodeParts[i]=NUM2INT(rb_ary_entry(rb_nodeParts,i));
  }
  faceId = NUM2INT(rb_faceId);
  for (i=0;i<6;i++) uvs[i] = NUM2DBL(rb_ary_entry(rb_uvs,i));
  return (queue==queueAddFace(queue,nodes,faceId,nodeParts,uvs)?self:Qnil);
}

VALUE queue_addedFaces( VALUE self, VALUE transaction )
{
  GET_QUEUE_FROM_SELF;
  return INT2NUM(queueAddedFaces(queue,NUM2INT(transaction)));
}

VALUE queue_addedFaceNodes( VALUE self, VALUE index )
{
  int i, nodes[3];
  VALUE rb_nodes;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedFaceNodes(queue,NUM2INT(index),nodes)) return Qnil;
  rb_nodes = rb_ary_new2(3);
  for (i=0;i<3;i++) rb_ary_store(rb_nodes,i,INT2NUM(nodes[i]));
  return rb_nodes;
}

VALUE queue_addedFaceId( VALUE self, VALUE index )
{
  int faceId;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedFaceId(queue,NUM2INT(index),&faceId)) return Qnil;
  return INT2NUM(faceId);
}

VALUE queue_addedFaceNodeParts( VALUE self, VALUE index )
{
  int i, nodeParts[3];
  VALUE rb_nodeParts;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedFaceNodeParts(queue,NUM2INT(index),nodeParts)) return Qnil;
  rb_nodeParts = rb_ary_new2(3);
  for (i=0;i<3;i++) rb_ary_store(rb_nodeParts,i,INT2NUM(nodeParts[i]));
  return rb_nodeParts;
}

VALUE queue_addedFaceUVs( VALUE self, VALUE index )
{
  int i;
  double uvs[6];
  VALUE rb_uvs;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedFaceUVs(queue,NUM2INT(index),uvs)) return Qnil;
  rb_uvs = rb_ary_new2(6);
  for (i=0;i<6;i++) rb_ary_store(rb_uvs,i,rb_float_new(uvs[i]));
  return rb_uvs;
}

/* ****************************** edges ****************************** */

VALUE queue_removeEdge( VALUE self, VALUE rb_nodes, VALUE rb_nodeParts )
{
  int i, nodes[2], nodeParts[2];
  GET_QUEUE_FROM_SELF;
  for (i=0;i<2;i++) nodes[i]=NUM2INT(rb_ary_entry(rb_nodes,i));
  for (i=0;i<2;i++) nodeParts[i]=NUM2INT(rb_ary_entry(rb_nodeParts,i));
  return (queue==queueRemoveEdge(queue,nodes,nodeParts)?self:Qnil);
}

VALUE queue_removedEdges( VALUE self, VALUE transaction )
{
  GET_QUEUE_FROM_SELF;
  return INT2NUM(queueRemovedEdges(queue,NUM2INT(transaction)));
}

VALUE queue_removedEdgeNodes( VALUE self, VALUE index )
{
  int i, nodes[2];
  VALUE rb_nodes;
  GET_QUEUE_FROM_SELF;
  if (queue != queueRemovedEdgeNodes(queue,NUM2INT(index),nodes)) return Qnil;
  rb_nodes = rb_ary_new2(2);
  for (i=0;i<2;i++) rb_ary_store(rb_nodes,i,INT2NUM(nodes[i]));
  return rb_nodes;
}

VALUE queue_removedEdgeNodeParts( VALUE self, VALUE index )
{
  int i, nodeParts[2];
  VALUE rb_nodeParts;
  GET_QUEUE_FROM_SELF;
  if (queue != queueRemovedEdgeNodeParts(queue,NUM2INT(index),nodeParts)) 
    return Qnil;
  rb_nodeParts = rb_ary_new2(2);
  for (i=0;i<2;i++) rb_ary_store(rb_nodeParts,i,INT2NUM(nodeParts[i]));
  return rb_nodeParts;
}

VALUE queue_addEdge( VALUE self, VALUE rb_nodes, VALUE rb_edgeId, 
		     VALUE rb_nodeParts, VALUE rb_ts )
{
  int i, nodes[2], edgeId, nodeParts[2];
  double ts[2];
  GET_QUEUE_FROM_SELF;
  for (i=0;i< 2;i++) {
    nodes[i]=NUM2INT(rb_ary_entry(rb_nodes,i));
    nodeParts[i]=NUM2INT(rb_ary_entry(rb_nodeParts,i));
  }
  edgeId = NUM2INT(rb_edgeId);
  for (i=0;i<2;i++) ts[i] = NUM2DBL(rb_ary_entry(rb_ts,i));
  return (queue==queueAddEdge(queue,nodes,edgeId,nodeParts,ts)?self:Qnil);
}

VALUE queue_addedEdges( VALUE self, VALUE transaction )
{
  GET_QUEUE_FROM_SELF;
  return INT2NUM(queueAddedEdges(queue,NUM2INT(transaction)));
}

VALUE queue_addedEdgeNodes( VALUE self, VALUE index )
{
  int i, nodes[2];
  VALUE rb_nodes;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedEdgeNodes(queue,NUM2INT(index),nodes)) return Qnil;
  rb_nodes = rb_ary_new2(2);
  for (i=0;i<2;i++) rb_ary_store(rb_nodes,i,INT2NUM(nodes[i]));
  return rb_nodes;
}

VALUE queue_addedEdgeId( VALUE self, VALUE index )
{
  int edgeId;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedEdgeId(queue,NUM2INT(index),&edgeId)) return Qnil;
  return INT2NUM(edgeId);
}

VALUE queue_addedEdgeNodeParts( VALUE self, VALUE index )
{
  int i, nodeParts[2];
  VALUE rb_nodeParts;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedEdgeNodeParts(queue,NUM2INT(index),nodeParts)) return Qnil;
  rb_nodeParts = rb_ary_new2(2);
  for (i=0;i<2;i++) rb_ary_store(rb_nodeParts,i,INT2NUM(nodeParts[i]));
  return rb_nodeParts;
}

VALUE queue_addedEdgeTs( VALUE self, VALUE index )
{
  int i;
  double ts[2];
  VALUE rb_ts;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedEdgeTs(queue,NUM2INT(index),ts)) return Qnil;
  rb_ts = rb_ary_new2(2);
  for (i=0;i<2;i++) rb_ary_store(rb_ts,i,rb_float_new(ts[i]));
  return rb_ts;
}

VALUE queue_dump( VALUE self )
{
  VALUE array;
  int i, nInt, nDouble;
  int *ints;
  double *doubles;
  GET_QUEUE_FROM_SELF;

  if (queue != queueDumpSize(queue,&nInt,&nDouble)) return Qnil;

  ints = malloc(nInt * sizeof(int));
  doubles = malloc(nDouble * sizeof(double));

  if (queue != queueDump(queue,ints,doubles)) {
    free(ints);
    free(doubles);
    return Qnil;
  }

  array = rb_ary_new2(nInt+nDouble);
  for (i=0;i<nInt;i++) rb_ary_store(array,i,INT2NUM(ints[i]));
  for (i=0;i<nDouble;i++) rb_ary_store(array,nInt+i,rb_float_new(doubles[i]));
  return array;
}

VALUE queue_dumpInt( VALUE self )
{
  VALUE array;
  int i, nInt, nDouble;
  int *ints;
  double *doubles;
  GET_QUEUE_FROM_SELF;

  if (queue != queueDumpSize(queue,&nInt,&nDouble)) return Qnil;

  ints = malloc(nInt * sizeof(int));
  doubles = malloc(nDouble * sizeof(double));

  if (queue != queueDump(queue,ints,doubles)) {
    free(ints);
    free(doubles);
    return Qnil;
  }

  array = rb_ary_new2(nInt);
  for (i=0;i<nInt;i++) rb_ary_store(array,i,INT2NUM(ints[i]));
  return array;
}

VALUE queue_dumpFloat( VALUE self )
{
  VALUE array;
  int i, nInt, nDouble;
  int *ints;
  double *doubles;
  GET_QUEUE_FROM_SELF;

  if (queue != queueDumpSize(queue,&nInt,&nDouble)) return Qnil;

  ints = malloc(nInt * sizeof(int));
  doubles = malloc(nDouble * sizeof(double));

  if (queue != queueDump(queue,ints,doubles)) {
    free(ints);
    free(doubles);
    return Qnil;
  }

  array = rb_ary_new2(nDouble);
  for (i=0;i<nDouble;i++) rb_ary_store(array,i,rb_float_new(doubles[i]));
  return array;
}

VALUE queue_load( VALUE self, VALUE rb_ints, VALUE rb_doubles )
{
  Queue *result;
  int i, nInt, nDouble;
  int *ints;
  double *doubles;
  GET_QUEUE_FROM_SELF;

  nInt = RARRAY(rb_ints)->len;
  nDouble = RARRAY(rb_doubles)->len;

  ints = malloc(nInt * sizeof(int));
  doubles = malloc(nDouble * sizeof(double));

  for (i=0;i<nInt;i++) ints[i] = NUM2INT(rb_ary_entry(rb_ints,i));
  for (i=0;i<nDouble;i++) doubles[i] = NUM2DBL(rb_ary_entry(rb_doubles,i));

  result = queueLoad( queue, ints, doubles );

  free(ints);
  free(doubles);

  return (queue==result?self:Qnil);
}

VALUE cQueue;

void Init_Queue() 
{
  cQueue = rb_define_class( "Queue", rb_cObject );
  rb_define_singleton_method( cQueue, "new", queue_new, 1 );
  rb_define_method( cQueue, "reset", queue_reset, 0 );
  rb_define_method( cQueue, "nodeSize", queue_nodeSize, 0 );
  rb_define_method( cQueue, "transactions", queue_transactions, 0 );
  rb_define_method( cQueue, "newTransaction", queue_newTransaction, 0 );

  rb_define_method( cQueue, "removeCell", queue_removeCell, 2 );
  rb_define_method( cQueue, "removedCells", queue_removedCells, 1 );
  rb_define_method( cQueue, "removedCellNodes", queue_removedCellNodes, 1 );
  rb_define_method( cQueue, "removedCellNodeParts", queue_removedCellNodeParts, 1 );
  rb_define_method( cQueue, "addCell", queue_addCell, 4 );
  rb_define_method( cQueue, "addedCells", queue_addedCells, 1 );
  rb_define_method( cQueue, "addedCellNodes", queue_addedCellNodes, 1 );
  rb_define_method( cQueue, "addedCellId", queue_addedCellId, 1 );
  rb_define_method( cQueue, "addedCellNodeParts", queue_addedCellNodeParts, 1 );
  rb_define_method( cQueue, "addedCellXYZs", queue_addedCellXYZs, 1 );
  rb_define_method( cQueue, "totalRemovedCells", queue_totalRemovedCells, 0 );
  
  rb_define_method( cQueue, "removeFace", queue_removeFace, 2 );
  rb_define_method( cQueue, "removedFaces", queue_removedFaces, 1 );
  rb_define_method( cQueue, "removedFaceNodes", queue_removedFaceNodes, 1 );
  rb_define_method( cQueue, "removedFaceNodeParts", queue_removedFaceNodeParts, 1 );
  rb_define_method( cQueue, "addFace", queue_addFace, 4 );
  rb_define_method( cQueue, "addedFaces", queue_addedFaces, 1 );
  rb_define_method( cQueue, "addedFaceNodes", queue_addedFaceNodes, 1 );
  rb_define_method( cQueue, "addedFaceId", queue_addedFaceId, 1 );
  rb_define_method( cQueue, "addedFaceNodeParts", queue_addedFaceNodeParts, 1 );
  rb_define_method( cQueue, "addedFaceUVs", queue_addedFaceUVs, 1 );

  rb_define_method( cQueue, "removeEdge", queue_removeEdge, 2 );
  rb_define_method( cQueue, "removedEdges", queue_removedEdges, 1 );
  rb_define_method( cQueue, "removedEdgeNodes", queue_removedEdgeNodes, 1 );
  rb_define_method( cQueue, "removedEdgeNodeParts", queue_removedEdgeNodeParts, 1 );
  rb_define_method( cQueue, "addEdge", queue_addEdge, 4 );
  rb_define_method( cQueue, "addedEdges", queue_addedEdges, 1 );
  rb_define_method( cQueue, "addedEdgeNodes", queue_addedEdgeNodes, 1 );
  rb_define_method( cQueue, "addedEdgeId", queue_addedEdgeId, 1 );
  rb_define_method( cQueue, "addedEdgeNodeParts", queue_addedEdgeNodeParts, 1 );
  rb_define_method( cQueue, "addedEdgeTs", queue_addedEdgeTs, 1 );

  rb_define_method( cQueue, "dump", queue_dump, 0 );  
  rb_define_method( cQueue, "dumpInt", queue_dumpInt, 0 );  
  rb_define_method( cQueue, "dumpFloat", queue_dumpFloat, 0 );  
  rb_define_method( cQueue, "load", queue_load, 2 );
}
