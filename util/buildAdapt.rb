#!/usr/bin/env ruby

class Package

 attr :path

 def initialize name
  @name = name
  @repository = case @name
		when "CAPRI", "SDK"
		 "-d :ext:geolab:/usr/local/inhouse/CVS"
		when "refine"
		 "-d :ext:cmb20:/ump/fldmd/home/mikepark/cvsroot"
		when "HEFSS.rps"
		 "-d :ext:hefss-core:/usr/local/cvsroot"
		else
		 puts "repository unknown: "+name
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

 def libs
  case @name
  when "SDK"
   "-L#{File.join(@path,'lib')} -lCADGeom-CAPRI -lMeat"
  when "CAPRI"
   "-L#{File.join(@path,'Native','LINUX')} -lcapri -L/usr/X11R6/lib -lX11"
  when "refine"
   "-L#{File.join(@path,'lib')} -lrefine"
  else
   puts "repository unknown: "+name
   nil
  end
 end

 def configure_env(capri_libs, sdk_libs, refine_libs)
  `(cd #@path && ./Configure)`
  File.open(File.join(@path,'Makefile.env'),'a') do |f|
   f.puts "CAPRILIBS = #{sdk_libs} #{capri_libs}"
   f.puts "REFINELIBS = #{refine_libs}"
   f.puts "F90FLAGS = -O0"
  end
  self
 end

 def make_make_install
  `(cd #@path && make > build.make && make install > build.install )`
  self
 end

 def make_mpi
  `(cd #@path && make mpi > build.make )`
  self
 end

 def adapt hefss_path
  `(cd #{File.join(@path,'example')} && ./RunAdapt.rb #{hefss_path} )`
 end

end

capri = Package.new("CAPRI")
sdk = Package.new("SDK")
refine = Package.new("refine")
hefss = Package.new("HEFSS.rps")

if false
 capri.clean_checkout
 sdk.clean_checkout.bootstrap.configure(capri.with).make_make_install
 refine.clean_checkout.bootstrap.configure(sdk.with,capri.with).make_make_install
 hefss.clean_checkout.configure_env(capri.libs,sdk.libs,refine.libs).make_mpi
end

refine.adapt(hefss.path)
