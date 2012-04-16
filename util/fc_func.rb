#!/usr/bin/env ruby

lines = IO.readlines(ARGV[0])

File.open(ARGV[0],'wt') do |f|
  lines.each do |line|
    if line =~ /void (.*)\(/
      func_ = $1
      func = $1.gsub(/_$/,'')
      puts line.sub!(/#{func_}/,"FC_FUNC_(#{func},#{func.upcase})")
    end
    f.puts line
  end
end
