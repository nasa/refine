#!/usr/bin/env ruby

starting_path = Dir.pwd

refine_path = "#{File.dirname $0}/../src"

$:.push refine_path

Dir.chdir refine_path
require 'RubyExtensionBuilder'

RubyExtensionBuilder.new('Grid').build
require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'

Dir.chdir starting_path

filename = ARGV[0]
grid = Grid.from_FAST File.expand_path(filename)
grid.exportGRI
grid.writeTecplotSurfaceGeom
