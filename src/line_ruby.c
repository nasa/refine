
#include "ruby.h"
#include "line.h"

#define GET_LINE_FROM_SELF Line *line; Data_Get_Struct( self, Line, line );

static void line_free( void *line )
{
  lineFree( line );
}

VALUE line_new( VALUE class )
{
  Line *line;
  VALUE obj;
  line = lineCreate( );
  obj = Data_Wrap_Struct( class, 0, line_free, line );
  return obj;
}

VALUE line_length( VALUE self )
{
  GET_LINE_FROM_SELF;
  return INT2NUM(lineLength(line));
}

VALUE line_addNode( VALUE self, VALUE node )
{
  GET_LINE_FROM_SELF;
  return (line==lineAddNode(line,NUM2INT(node))?self:Qnil);
}

VALUE line_node( VALUE self, VALUE index )
{
  GET_LINE_FROM_SELF;
  return INT2NUM(lineNode(line,NUM2INT(index)));
}

VALUE cLine;

#define GET_LINES_FROM_SELF Lines *lines; Data_Get_Struct( self, Lines, lines);

static void lines_free( void *lines )
{
  linesFree( lines );
}

VALUE lines_new( VALUE class )
{
  Lines *lines;
  VALUE obj;
  lines = linesCreate( );
  obj = Data_Wrap_Struct( class, 0, lines_free, lines );
  return obj;
}

VALUE lines_number( VALUE self )
{
  GET_LINES_FROM_SELF;
  return INT2NUM(linesNumber(lines));
}

VALUE lines_addNode( VALUE self, VALUE line, VALUE node )
{
  GET_LINES_FROM_SELF;
  return (lines==linesAddNode(lines,NUM2INT(line),NUM2INT(node))?self:Qnil);
}

VALUE lines_node( VALUE self, VALUE line, VALUE index )
{
  GET_LINES_FROM_SELF;
  return INT2NUM(linesNode(lines,NUM2INT(line),NUM2INT(index)));
}

VALUE lines_renumber( VALUE self, VALUE rb_o2n )
{
  int i, length, *o2n;
  GET_LINES_FROM_SELF;
  
  length = RARRAY_LEN(rb_o2n);
  o2n = malloc( length * sizeof(int) );
  for (i=0;i<length;i++) o2n[i] = NUM2INT(rb_ary_entry(rb_o2n,i));
  linesRenumber(lines,o2n);
  free(o2n);
  return self;
}

VALUE lines_save( VALUE self, VALUE rb_filename )
{
  char *filename;
  GET_LINES_FROM_SELF;
  filename = RSTRING(rb_filename)->ptr;
  return (lines==linesSave(lines,filename)?self:Qnil);
}

VALUE lines_load( VALUE self, VALUE rb_filename )
{
  char *filename;
  GET_LINES_FROM_SELF;
  filename = RSTRING(rb_filename)->ptr;
  return (lines==linesLoad(lines,filename)?self:Qnil);
}

VALUE cLines;

void Init_Line() 
{
  cLine = rb_define_class( "Line", rb_cObject );
  rb_define_singleton_method( cLine, "new", line_new, 0 );
  rb_define_method( cLine, "length", line_length, 0 );
  rb_define_method( cLine, "addNode", line_addNode, 1 );
  rb_define_method( cLine, "node", line_node, 1 );

  cLines = rb_define_class( "Lines", rb_cObject );
  rb_define_singleton_method( cLines, "new", lines_new, 0 );
  rb_define_method( cLines, "number", lines_number, 0 );
  rb_define_method( cLines, "addNode", lines_addNode, 2 );
  rb_define_method( cLines, "node", lines_node, 2 );
  rb_define_method( cLines, "renumber", lines_renumber, 1 );
  rb_define_method( cLines, "save", lines_save, 1 );
  rb_define_method( cLines, "load", lines_load, 1 );
}
