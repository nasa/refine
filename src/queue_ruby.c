
#include "ruby.h"
#include "queue.h"

#define GET_QUEUE_FROM_SELF Queue *queue; Data_Get_Struct( self, Queue, queue );

static void queue_free( void *queue )
{
  queueFree( queue );
}

VALUE queue_new( VALUE class )
{
  Queue *queue;
  VALUE obj;
  queue = queueCreate(  );
  obj = Data_Wrap_Struct( class, 0, queue_free, queue );
  return obj;
}

VALUE queue_reset(VALUE self )
{
  GET_QUEUE_FROM_SELF;
  return (queue==queueReset(queue)?self:Qnil);
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

VALUE queue_removeCell( VALUE self, VALUE rb_nodes )
{
  int i, nodes[4];
  GET_QUEUE_FROM_SELF;
  for (i=0;i<4;i++) nodes[i]=NUM2INT(rb_ary_entry(rb_nodes,i));
  return (queue==queueRemoveCell(queue,nodes)?self:Qnil);
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

VALUE queue_addCell( VALUE self, VALUE rb_nodes, VALUE rb_xyzs )
{
  int i, nodes[5];
  double xyzs[12];
  GET_QUEUE_FROM_SELF;
  for (i=0;i< 5;i++) nodes[i]=NUM2INT(rb_ary_entry(rb_nodes,i));
  for (i=0;i<12;i++) xyzs[i] =NUM2DBL(rb_ary_entry(rb_xyzs,i));
  return (queue==queueAddCell(queue,nodes,xyzs)?self:Qnil);
}

VALUE queue_addedCells( VALUE self, VALUE transaction )
{
  GET_QUEUE_FROM_SELF;
  return INT2NUM(queueAddedCells(queue,NUM2INT(transaction)));
}

VALUE queue_addedCellNodes( VALUE self, VALUE index )
{
  int i, nodes[5];
  VALUE rb_nodes;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedCellNodes(queue,NUM2INT(index),nodes)) return Qnil;
  rb_nodes = rb_ary_new2(5);
  for (i=0;i<5;i++) rb_ary_store(rb_nodes,i,INT2NUM(nodes[i]));
  return rb_nodes;
}

VALUE queue_addedCellXYZs( VALUE self, VALUE index )
{
  int i;
  double xyzs[12];
  VALUE rb_xyzs;
  GET_QUEUE_FROM_SELF;
  if (queue != queueAddedCellXYZs(queue,NUM2INT(index),xyzs)) return Qnil;
  rb_xyzs = rb_ary_new2(12);
  for (i=0;i<12;i++) rb_ary_store(rb_xyzs,i,rb_float_new(xyzs[i]));
  return rb_xyzs;
}

VALUE cQueue;

void Init_Queue() 
{
  cQueue = rb_define_class( "Queue", rb_cObject );
  rb_define_singleton_method( cQueue, "new", queue_new, 0 );
  rb_define_method( cQueue, "reset", queue_reset, 0 );
  rb_define_method( cQueue, "transactions", queue_transactions, 0 );
  rb_define_method( cQueue, "newTransaction", queue_newTransaction, 0 );
  rb_define_method( cQueue, "removeCell", queue_removeCell, 1 );
  rb_define_method( cQueue, "removedCells", queue_removedCells, 1 );
  rb_define_method( cQueue, "removedCellNodes", queue_removedCellNodes, 1 );
  rb_define_method( cQueue, "addCell", queue_addCell, 2 );
  rb_define_method( cQueue, "addedCells", queue_addedCells, 1 );
  rb_define_method( cQueue, "addedCellNodes", queue_addedCellNodes, 1 );
  rb_define_method( cQueue, "addedCellXYZs", queue_addedCellXYZs, 1 );
}
