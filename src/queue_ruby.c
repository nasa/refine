
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

VALUE queue_transactionNodes( VALUE self, VALUE transaction )
{
  GET_QUEUE_FROM_SELF;
  return INT2NUM(queueTransactionNodes(queue,NUM2INT(transaction)));
}

VALUE queue_addNode( VALUE self, VALUE node )
{
  GET_QUEUE_FROM_SELF;
  return (queue==queueAddNode(queue,NUM2INT(node))?self:Qnil);
}

VALUE cQueue;

void Init_Queue() 
{
  cQueue = rb_define_class( "Queue", rb_cObject );
  rb_define_singleton_method( cQueue, "new", queue_new, 0 );
  rb_define_method( cQueue, "reset", queue_reset, 0 );
  rb_define_method( cQueue, "transactions", queue_transactions, 0 );
  rb_define_method( cQueue, "newTransaction", queue_newTransaction, 0 );
  rb_define_method( cQueue, "transactionNodes", queue_transactionNodes, 1 );
  rb_define_method( cQueue, "addNode", queue_addNode, 1 );
}
