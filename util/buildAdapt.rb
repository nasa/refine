#!/usr/bin/env ruby

class Package

 def initialize name
  @name = name
  @repository = case name
		when "CAPRI", "SDK"
		 "-d :ext:geolab:/usr/local/inhouse/CVS"
		when "refine"
		 "-d :ext:cmb20:/ump/fldmd/home/mikepark/cvsroot"
		when "HEFSS.rps"
		 "-d :ext:hefss-core:/usr/local/cvsroot"
		else
		 puts "repository unknown: "+package
		 nil
		end
  @where=`pwd`
 end

 def kill
  `rm -rf #@name`
  self
 end

 def check_out
  command = "cvs -q #@repository co -P #@name"
  puts command
  `#{command}`
  self
 end

 def clean_check_out
  kill
  check_out
  self
 end
end

capri = Package.new("CAPRI")
sdk = Package.new("SDK")
refine = Package.new("refine")
hefss = Package.new("HEFSS.rps")
