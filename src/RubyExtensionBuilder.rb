# Ruby C extension build for refine package
#
# $Id$

class RubyExtensionBuilder

 def initialize extension
  @extension = extension
 end

 def buildOnly extension
  comm = Hash.new('puts "Error, extension not known";exit 1')
  comm['Adj'] = 'Adj'
  comm['Near'] = 'Near'
  comm['Intersect'] = 'Intersect gridmath.c'
  comm['Grid'] = 'Grid adj.h'
  comm['GridMath'] = 'GridMath'
  comm['GridMetric'] = 'GridMetric adj.h grid.h gridmath.c'
  comm['GridSwap'] = 'GridSwap adj.h grid.h gridmath.h gridmetric.h'
  comm['GridCAD'] = 'GridCAD FAKEGeom adj.h grid.h gridmath.h gridmetric.h gridinsert.h'
  comm['GridInsert'] = 'GridInsert adj.h grid.h gridmath.h gridmetric.h gridcad.h'
  comm['Layer'] = 'Layer layerStruct.h adj.h near.h intersect.h grid.h gridmath.h gridmetric.h gridcad.h gridinsert.h'
  exit 1 unless system('ruby makeRubyExtension.rb ' + comm[extension] + ' master_header.h')
 end

 def build
  depend = Hash.new
  depend['Adj'] = %w{ Adj }
  depend['Near'] = %w{ Near }
  depend['Intersect'] = %w{ Intersect }
  depend['Grid'] = %w{ Adj Grid }
  depend['GridMath'] = %w{ GridMath }
  depend['GridMetric'] = %w{ Adj Grid GridMetric }
  depend['GridSwap'] = %w{ Adj Grid GridMetric GridSwap }
  depend['GridCAD'] = %w{ Adj Grid GridMetric GridCAD }
  depend['GridInsert'] = %w{ Adj Grid GridMetric GridSwap GridCAD GridInsert }
  depend['Layer'] = %w{ Adj Near Grid GridMetric GridCAD GridInsert Layer }
  depend[@extension].each { |extension|
   buildOnly extension
  }
 end

end
