#  node, x-y pair

class Node
 
 attr_accessor :x, :y, :bc, :i

 def initialize(x, y)
  self.x = x
  self.y = y
 end

 def to_s
   "#{x} #{y} #{bc}\n"
 end

end
#  

class Segment
 
 attr_accessor :bc, :first, :last

 def initialize(bc,first,last)
  self.bc = bc
  self.first = first
  self.last = last
  first.bc = bc
 end

 def to_s
   "#{first.i} #{last.i} #{bc}\n"
 end

end
#  node, x-y pair

class Msh < Array

 def header
  "Dimension\n2\nMaximalAngleOfCorner\n0.001\n"
 end

 def set_node_i_based_on_segment_list
  i = 0
  self.each do |segment|
   node = segment.first
   i += 1
   node.i = i
  end
 end

 def vertices
  set_node_i_based_on_segment_list
  n = self.length
  out = "Vertices\n#{n}\n"
  self.each do |segment|
   node = segment.first
   out += node.to_s
  end
  out
 end

 def edges
  n = self.length
  out = "Edges\n#{n}\n"
  self.each do |segment|
   out += segment.to_s
  end
  out
 end

 def subdomain
   "SubDomain\n1\n2 1 1 0"
 end

 def to_s
  header+vertices+edges+subdomain
 end

 def export(filename='msh.msh')
  File.open(filename,'w') do |f|
   f.puts to_s
  end
 end
end
