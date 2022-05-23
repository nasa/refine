#!/usr/bin/env ruby

def decore(filename)
puts filename
lines = IO.readlines(filename)

File.open(filename,'w') do |f|
  lines.each do |line|
    line.chomp!
    line.sub!(/^REF_STATUS ref/, "REF_FCN REF_STATUS ref")
    line.sub!(/^static REF_STATUS ref/, "REF_FCN static REF_STATUS ref")
    line.sub!(/^REF_BOOL ref/, "REF_FCN REF_BOOL ref")
    line.sub!(/^static REF_BOOL ref/, "REF_FCN static REF_BOOL ref")
    f.puts line
  end
end
end

ARGV.each do |filename|
  decore(filename)
end


