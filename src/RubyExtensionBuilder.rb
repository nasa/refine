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
  extraFiles['Grid'] = 'adj.h line.h queue.h sort.h'
  extraFiles['GridMetric'] = 'adj.h line.h queue.h grid.h gridmath.h'
  extraFiles['GridSwap'] = 'adj.h line.h queue.h grid.h gridmath.h gridmetric.h'
  extraFiles['GridCAD'] = 'FAKEGeom.c adj.h line.h grid.h gridmath.h gridmetric.h gridswap.h queue.h gridinsert.h'
  extraFiles['GridInsert'] = 'adj.h line.h grid.h gridmath.h gridswap.h gridmetric.h gridcad.h queue.h '
  extraFiles['GridMPI'] = 'adj.h line.h queue.h grid.h gridmath.h gridmetric.h gridinsert.h gridswap.h'
  extraFiles['Layer'] = 'layerStruct.h adj.h line.h grid.h gridmath.h near.h intersect.h gridmetric.h gridcad.h queue.h gridinsert.h'

  systemCall = ['ruby makeRubyExtension.rb',extension,extraFiles[extension],'refine_defs.h'].join(' ')
  exit 1 unless system(systemCall )
 end

 def build
  requiredPackages = Hash.new([])
  requiredPackages['Grid'] = %w{ Adj Line Sort}
  requiredPackages['GridMetric'] = %w{ Adj Line Grid }
  requiredPackages['GridSwap'] = %w{ Adj Line Grid GridMath GridMetric }
  requiredPackages['GridCAD'] = %w{ Adj Line Grid GridMath GridMetric }
  requiredPackages['GridInsert'] = %w{ Adj Line Grid GridMetric GridSwap GridCAD }
  requiredPackages['GridMPI'] = %w{ Adj Line Sort Queue Grid GridMath GridMetric GridInsert GridSwap }
  requiredPackages['Layer'] = %w{ Adj Near Intersect Line Grid GridMath GridMetric GridCAD GridInsert }

  requiredPackages[@extension].each { |extension| buildOnly extension }
  buildOnly @extension
 end

end
