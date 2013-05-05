
require 'scanf'

#  node, x-y pair

class Node
 
 attr_accessor :x, :y, :bc, :i

 def initialize(x, y, bc = nil)
  self.x = x
  self.y = y
  self.bc = bc
 end

 def to_s
   "#{x} #{y} #{bc}\n"
 end

end


class Curve < Array

  attr_accessor :bc

  def initialize(bc,filename, scale=1.0, format="%f,%f")
    self.bc = bc

    File.open(filename).each do |line|
      values = line.scanf(format)
      if (values.length >= 2)
        if ( !last.nil? &&
             values[0]*scale == last.x &&
             values[1]*scale == last.y ) 
          puts "duplicate node dropped"
        else
          self << Node.new(values[0]*scale,values[1]*scale,bc) 
        end
      end
    end

  end

  def to_s
    out=""
    (self.length-1).times do |i|
      out += "#{self[i].i} #{self[i+1].i} #{bc}\n"
    end
    out
  end

end


class Segment < Array
 
 attr_accessor :bc

 def initialize(bc,node0,node1)
  self.bc = bc
  self[0] = node0
  self[1] = node1
  first.bc = bc
 end

 def to_s
   "#{first.i} #{last.i} #{bc}\n"
 end

end


class Msh < Array

 def initialize(angle=30)
   @verts = Array.new
   @angle=angle
 end

 def header
   "AngleOfCornerBound\n#{@angle}\nDimension\n2\n"
 end

 def set_node_i_based_on_segment_list
   i = 0
   self.each do |segment|
     segment.each do |node|
       unless node.i
         i += 1
         node.i = i
         @verts << node
       end
     end
   end
 end

 def vertices
  set_node_i_based_on_segment_list
  n = @verts.length
  out = "Vertices\n#{n}\n"
  @verts.each do |node|
   out += node.to_s
  end
  out
 end

 def edges
   n = 0
   self.each do |segment|
     n += (segment.length-1)
   end
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
