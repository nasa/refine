
# do all if no args and put in rubyext

require 'mkmf'


unless ext = ARGV[0] 
  puts "ERROR: usage: ruby extconf.rb rubyExtensionName"
  exit 1
end



rubyExt = ext + "_ruby.c"

sourceExt = rubyExt.collect{ |r| r.sub(/_ruby.c/, ".c") }

objC = rubyExt.to_a + sourceExt

$objs = objC.collect{ |c| c.sub(/.c/, ".o") }

create_makefile(ext)
