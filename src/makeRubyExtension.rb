#!/usr/bin/env ruby
#
# Creates a ruby c extension from the args
#
# $Id$

unless ext = ARGV[0] 
 puts "ERROR: usage: ruby #{__FILE__} rubyExtensionName [extraFiles.(c|h)]"
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
 require 'mkmf'
 $objs = objC.collect{ |c| c.sub(/\.c/, ".o") }
 create_makefile(ext,'..')
 exit 1 unless system "make --quiet --no-print-directory"
Dir.chdir '..'
