#!/usr/bin/env ruby

files = Dir["*.h"]+Dir["*.c"]

files.each do |file|
  contents = IO.readlines(file)
  File.open(file,'w') do |f|
    contents.each do |line|
      f.puts line.gsub(/#{ARGV[0]}/,ARGV[1])
    end
  end
end
