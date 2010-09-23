#!/usr/bin/env ruby

refine_path = File.expand_path("~/GRIDEX/refine/src")

$:.push refine_path

require 'Adj/Adj'
require 'Line/Line'
require 'Grid/Grid'

Grid.from_FAST(ARGV[0]).exportFASTSurface

