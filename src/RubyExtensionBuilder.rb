# Ruby C extension build for refine package
#
# $Id$

class RubyExtensionBuilder

 def initialize extension
  @extension = extension
 end

 def buildOnly extension
  extraFiles = Hash.new('')
  extraFiles['Intersect'] = 'gridmath.c'
  extraFiles['Grid'] = 'adj.h'
  extraFiles['GridMetric'] = 'adj.h grid.h gridmath.c'
  extraFiles['GridSwap'] = 'adj.h grid.h gridmath.h gridmetric.h'
  extraFiles['GridCAD'] = 'FAKEGeom adj.h grid.h gridmath.h gridmetric.h gridinsert.h'
  extraFiles['GridInsert'] = 'adj.h grid.h gridmath.h gridmetric.h gridcad.h'
  extraFiles['Layer'] = 'layerStruct.h adj.h near.h intersect.h grid.h gridmath.h gridmetric.h gridcad.h gridinsert.h'

  systemCall = ['ruby makeRubyExtension.rb',extension,extraFiles[extension],'master_header.h'].join(' ')
  exit 1 unless system(systemCall )
 end

 def build
  requiredPackages = Hash.new([])
  requiredPackages['Grid'] = %w{ Adj }
  requiredPackages['GridMetric'] = %w{ Adj Grid }
  requiredPackages['GridSwap'] = %w{ Adj Grid GridMetric }
  requiredPackages['GridCAD'] = %w{ Adj Grid GridMetric }
  requiredPackages['GridInsert'] = %w{ Adj Grid GridMetric GridSwap GridCAD }
  requiredPackages['Layer'] = %w{ Adj Near Intersect Grid GridMetric GridCAD GridInsert }

  requiredPackages[@extension].each { |extension| buildOnly extension }
  buildOnly @extension
 end

end
