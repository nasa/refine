
require 'mkmf'

rubyExt = Dir["*_ruby.c"]
sourceExt = rubyExt.collect{ |r| r.sub(/_ruby.c/, ".c") }


objC = rubyExt + sourceExt


$objs = objC.collect{ |c| c.sub(/.c/, ".o") }

create_makefile("SampleUnit")
