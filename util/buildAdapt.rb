#!/usr/bin/env ruby

def kill working_copy
 `rm -rf ${working_copy}`
end

def check_out package
 repository = case package
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
 `cvs -q ${repository} co -P ${package}`
end

