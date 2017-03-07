#! /usr/bin/env ruby

status_file = ARGV[0] || 'accept-3d-one-03.status' 
min_limit = ARGV[1] || '0.038'
max_limit = ARGV[2] || '2.3'
min_limit = min_limit.to_f
max_limit = max_limit.to_f

lines = IO.readlines(status_file)

min_edge = nil
max_edge = nil
lines.each do |line|
  min_edge = line.split.first.to_f if (line =~ /min edge/ )
  max_edge = line.split(':').first.to_f if (line =~ /max edge/ )
end

throw("min edge missing") if (min_edge.nil?)
throw("max edge missing") if (max_edge.nil?)

need_to_throw = false
if ( min_limit > min_edge )
  puts "min edge #{min_edge} short of limit #{min_limit}"
  need_to_throw = true
end
if ( max_limit < max_edge )
  puts "max edge #{max_edge} long of limit #{max_limit}"
  need_to_throw = true
end

throw("edge out of range") if (need_to_throw)

