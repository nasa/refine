
# do all if no args and put in rubyext

unless ext = ARGV[0] 
  puts "ERROR: usage: ruby extconf.rb rubyExtensionName"
  exit 1
end

rubyExt = ext.downcase + "_ruby.c"

sourceExt = rubyExt.collect{ |r| r.sub(/_ruby\.c/, ".c") }

objC = rubyExt.to_a + sourceExt

ARGV[1..ARGV.size].each { |c| objC.push c if c =~/\.c/ }

headers = objC.collect do
 |c| h = c.sub(/\.c/, ".h") 
 h if File.exist? h
end
ARGV[1..ARGV.size].each { |h| headers.push h if h =~/\.h/ }

`mkdir -p #{ext}`
Dir.chdir ext
(objC+headers).compact.each { |f| `ln -sf ../#{f} .`}

if ARGV.include?("FAKEGeom")
 `ln -sf ../FAKEGeom.c CADGeom.c`
 `mkdir -p CADGeom`
 Dir.chdir "CADGeom"
 `ln -sf ../../FAKEGeom.h CADGeom.h`
 `ln -sf ../../master_header.h .`
 Dir.chdir ".."
 objC.push "CADGeom.c"
end

$objs = objC.collect{ |c| c.sub(/\.c/, ".o") }

require 'mkmf'

create_makefile(ext)

exit 1 unless system "make"
