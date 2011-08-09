#!/usr/bin/env ruby

def extend_with( ext ) 

  rubyExt = ext.downcase + "_ruby.c"

  sourceExt = ext.downcase + ".c"

  objC = [ rubyExt, sourceExt ]

#  ARGV[1..ARGV.size].each { |c| objC.push c if c =~/\.c/ }

  headers = objC.collect do
    |c| h = c.sub(/\.c/, ".h") 
    h if File.exist? h
  end
#  ARGV[1..ARGV.size].each { |h| headers.push h if h =~/\.h/ }

  `chmod u+w . && mkdir -p #{ext}`
  Dir.chdir ext
  unless FileTest.exists?('Makefile')
    require 'mkmf'
    $objs = objC.collect{ |c| c.sub(/\.c/, ".o") }
    create_makefile(ext,'..')
  end
  exit 1 unless system "make --quiet --no-print-directory"
  Dir.chdir '..'

end

if (__FILE__ == $0) then 
  unless ext = ARGV[0] 
    puts "ERROR: usage: ruby #{__FILE__} rubyExtensionName [extraFiles.(c|h)]"
    exit 1
  end
  extend_with( ext )
end

