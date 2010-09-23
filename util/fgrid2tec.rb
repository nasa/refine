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

Grid.from_FAST(ARGV[0]).writeTecplotSurfaceGeom

