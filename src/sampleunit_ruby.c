
#include "ruby.h"
#include "sampleunit.h"

static VALUE su_sampleunit( VALUE self, VALUE f1, VALUE f2 )
{
  int sum;
  sum = sampleunit( FIX2INT( f1 ), FIX2INT( f2 ) );
  return INT2FIX( sum );
}

VALUE cSampleUnit;

void Init_SampleUnit(  )
{
  cSampleUnit = rb_define_class( "SampleUnit", rb_cObject );
  rb_define_method( cSampleUnit, "sampleUnit", su_sampleunit, 2 );
}
