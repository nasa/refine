
#include "ruby.h"
#include "sort.h"

static VALUE Sort_Heap( VALUE self, VALUE rb_arrayInput )
{
  int i, length;
  int *arrayInput, *sortedIndex;
  VALUE rb_sortedIndex;

  length = RARRAY(rb_arrayInput)->len;

  arrayInput = malloc(length*sizeof(int));
  sortedIndex = malloc(length*sizeof(int));
  for(i=0;i<length;i++) arrayInput[i]=NUM2INT(rb_ary_entry(rb_arrayInput,i));

  sortHeap( length, arrayInput, sortedIndex );

  rb_sortedIndex = rb_ary_new2(length);
  for(i=0;i<length;i++) rb_ary_store(rb_sortedIndex,i,INT2NUM(sortedIndex[i]));

  free(arrayInput);
  free(sortedIndex);

  return rb_sortedIndex;
}

static VALUE Sort_Search( VALUE self, VALUE rb_arrayInput, VALUE index )
{
  int i, length, found;
  int *arrayInput;

  length = RARRAY(rb_arrayInput)->len;

  arrayInput = malloc(length*sizeof(int));
  for(i=0;i<length;i++) arrayInput[i]=NUM2INT(rb_ary_entry(rb_arrayInput,i));

  found = sortSearch( length, arrayInput, index );

  free(arrayInput);

  return INT2NUM(found);
}

VALUE cSort;

void Init_Sort(  )
{
  cSort = rb_define_class( "Sort", rb_cObject );
  rb_define_singleton_method( cSort, "Heap", Sort_Heap, 1 );
  rb_define_singleton_method( cSort, "Search", Sort_Search, 2 );
}
