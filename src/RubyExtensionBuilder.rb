# Ruby C extension build for refine package
#
# $Id$

class RubyExtensionBuilder

 def initialize extension
  @extension = extension
 end

 def buildOnly extension
  extraFiles = Hash.new('')
  extraFiles['GridCAD'] = 'FAKEGeom.c'
  systemCall = ['ruby makeRubyExtension.rb',extension,extraFiles[extension],'refine_defs.h'].join(' ')
  exit 1 unless system(systemCall )
 end

 def build
  requiredPackages = Hash.new([])
  requiredPackages['Intersect'] = %w{ GridMath }
  requiredPackages['Grid'] = %w{ Adj Line Sort}
  requiredPackages['GridMetric'] = %w{ Adj Line Grid GridMath }
  requiredPackages['GridSwap'] = %w{ Adj Line Grid GridMath GridMetric }
  requiredPackages['GridCAD'] = %w{ Adj Line Grid GridMath GridMetric }
  requiredPackages['GridInsert'] = %w{ Adj Line Grid GridMetric GridSwap GridCAD }
  requiredPackages['GridMPI'] = %w{ Adj Line Sort Queue Grid GridMath GridMetric GridInsert GridSwap }
  requiredPackages['GridMove'] = %w{ Adj Line Grid GridMetric }
  requiredPackages['Layer'] = %w{ Adj Near Intersect Line Grid GridMath GridMetric GridCAD GridInsert }

  requiredPackages[@extension].each { |extension| buildOnly extension }
  buildOnly @extension
 end

end
