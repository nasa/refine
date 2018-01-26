#! /usr/bin/env ruby

status_file = ARGV[0] || 'recon' 
max_limit = ARGV[1] || '2e-12'
max_limit = max_limit.to_f

lines = IO.readlines(status_file)

max_node = nil
max_edge = nil
max_face = nil
lines.each do |line|
  max_node = line.split.last.to_f if (line =~ /node max eval dist/)
  max_edge = line.split.last.to_f if (line =~ /edge max eval dist/)
  max_face = line.split.last.to_f if (line =~ /face max eval dist/)
end

puts max_node
puts max_edge
puts max_face

throw("max node missing") if (max_node.nil?)
throw("max edge missing") if (max_edge.nil?)
throw("max face missing") if (max_face.nil?)

if ( max_limit < max_node )
  puts "max node #{max_node} outside of limit #{max_limit}"
  need_to_throw = true
end
if ( max_limit < max_edge )
  puts "max edge #{max_edge} outside of limit #{max_limit}"
  need_to_throw = true
end
if ( max_limit < max_face )
  puts "max face #{max_face} outside of limit #{max_limit}"
  need_to_throw = true
end

if (need_to_throw)
  puts `pwd`
  throw("projection out of range")
end
