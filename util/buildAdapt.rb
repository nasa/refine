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
  @where=`pwd`.chomp
  @path = File.join(@where,@name)
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

 def bootstrap
  autogen = File.join(path,"autogen.sh")
  `(cd #@path && #{autogen} > autogen.out)` if File.executable?(autogen)
 end

end

capri = Package.new("CAPRI")
sdk = Package.new("SDK")
refine = Package.new("refine")
hefss = Package.new("HEFSS.rps")

sdk.bootstrap
