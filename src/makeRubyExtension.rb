
# do all if no args and put in rubyext

unless ext = ARGV[0] 
  puts "ERROR: usage: ruby extconf.rb rubyExtensionName"
  exit 1
end

rubyExt = ext.downcase + "_ruby.c"

sourceExt = rubyExt.collect{ |r| r.sub(/_ruby.c/, ".c") }

objC = rubyExt.to_a + sourceExt
objC.push ARGV[1] if ARGV[1]

headers = objC.collect do
 |c| h = c.sub(/.c/, ".h") 
 h if File.exist? h
end
headers.push "master_header.h"

`rm -rf #{ext}`
`mkdir #{ext}`
`cp #{(objC+headers).join(' ')} #{ext}`
Dir.chdir ext

$objs = objC.collect{ |c| c.sub(/.c/, ".o") }

require 'mkmf'

create_makefile(ext)

system "make"
