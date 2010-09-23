# Ruby C extension build for refine package
#


class RubyExtensionBuilder

 def initialize extension
  @extension = extension
 end

 def buildOnly extension
  extraFiles = Hash.new('')
  extraFiles['GridCAD'] = 'FAKEGeom.c FAKEGeomExtras.c'
  extraFiles['GridMetric'] = 'FAKEGeom.c FAKEGeomExtras.c'
  systemCall = ['ruby makeRubyExtension.rb',extension,extraFiles[extension]].join(' ')
  exit 1 unless system( systemCall )
 end

 def build
  requiredPackages = Hash.new([])
  requiredPackages['Intersect'] = %w{ GridMath }
  requiredPackages['Plan'] = %w{ Sort }
  requiredPackages['Grid'] = %w{ Adj Line Sort GridMath }
  requiredPackages['GridShape'] = %w{ Adj Line Grid GridMath }
  requiredPackages['GridMetric'] = %w{ Adj Line Grid GridShape GridMath }
  requiredPackages['GridSwap'] = %w{ Adj Line Plan Sort Queue Grid GridMath GridShape GridMetric }
  requiredPackages['GridCAD'] = %w{ Adj Line Plan Sort Tableau Grid GridMath GridShape GridMetric }
  requiredPackages['GridInsert'] = %w{ Adj Line Intersect Plan Sort Queue Grid GridMath GridShape GridMetric GridSwap GridCAD }
  requiredPackages['GridMPI'] = %w{ Adj Line Sort Queue Grid GridMath GridShape GridMetric GridCAD GridInsert GridSwap }
  requiredPackages['GridMove'] = %w{ Adj Line Grid GridMath GridShape GridMetric }
  requiredPackages['GridEdger'] = %w{ Adj Line Queue Grid GridMath GridShape GridMetric GridCAD GridInsert}
  requiredPackages['GridFacer'] = %w{ Adj Line Grid GridMath GridShape GridMetric GridCAD GridInsert}
  requiredPackages['Layer'] = %w{ Adj Near Intersect Line Grid GridMath GridShape GridMetric GridCAD GridInsert }

  requiredPackages[@extension].each { |extension| buildOnly extension }
  buildOnly @extension
 end

end
