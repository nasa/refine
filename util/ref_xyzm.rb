#!/usr/bin/env ruby

ugrid = IO.readlines(ARGV[0])
metric = IO.readlines(ARGV[1])
output_filename = ARGV[2]

puts ugrid[0]
puts metric[0]

nnode = ugrid[0].split.first.to_i
puts nnode

File.open(output_filename,'w') do |f|
  nnode.times do |node|
    xyz = ugrid[node+1].chomp
    m = metric[node].chomp
    f.puts [xyz,m].join(' ')
  end
end
