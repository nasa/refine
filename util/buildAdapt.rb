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

 def checkout
  command = "cvs -q #@repository co -P #@name"
  puts command
  `#{command}`
  self
 end

 def clean_checkout
  kill
  checkout
  self
 end

 def bootstrap
  autogen = File.join(@path,"autogen.sh")
  `(cd #@path && #{autogen} > build.autogen)` if File.executable?(autogen)
  self
 end

 def with
  "--with-#{@name}=#{@path}"
 end

 def configure( option1="", option2="" ) 
  `(cd #@path && ./configure --prefix=#@path #{option1} #{option2} > build.configure)`
  self
 end

 def make_make_install
  `(cd #@path && make > build.make && make install > build.install )`
  self
 end

end

capri = Package.new("CAPRI")
sdk = Package.new("SDK")
refine = Package.new("refine")
hefss = Package.new("HEFSS.rps")

capri.clean_checkout
sdk.clean_checkout.bootstrap.configure(capri.with).make_make_install
refine.clean_checkout.bootstrap.configure(sdk.with,capri.with).make_make_install
