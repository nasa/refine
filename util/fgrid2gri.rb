#!/usr/bin/env ruby

refine_path = File.expand_path("~/GRIDEX/refine/src")

$:.push refine_path

# for Grid...
require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'

filename = ARGV[0]
grid = Grid.from_FAST File.expand_path(filename)
grid.exportGRI
grid.writeTecplotSurfaceGeom
