#!/usr/bin/env ruby
#
# Mobility test for line c lib
#


Dir.chdir ENV['srcdir'] if ENV['srcdir']
$:.push "." # ruby 1.9.2
require 'RubyExtensionBuilder'
RubyExtensionBuilder.new('Line').build

require 'test/unit'
require 'Line/Line'

class TestLine < Test::Unit::TestCase

 EMPTY = -1

 def testInitialLineLength
  line = Line.new
  assert_equal 0, line.length
 end

 def testInitializeEmptyNodes
  line = Line.new
  assert_equal EMPTY, line.node(-1)
  assert_equal EMPTY, line.node(0)
  assert_equal EMPTY, line.node(1)
 end

 def testAddOneNode
  line = Line.new
  node = 9
  line.addNode(node)
  assert_equal 1, line.length
  assert_equal node, line.node(0)
  assert_equal EMPTY, line.node(1)
 end

 def testAddSomeNodes
  line = Line.new
  nodes = [44, 47, 23, 87, 89, 3, 567]
  nodes.each { |node| line.addNode(node) }
  nodes.each_index do |index| 
   assert_equal nodes[index], line.node(index), "node #{index} doesnt match" 
  end
 end

 def testDontAddDuplicateNodes
  line = Line.new
  nodes = [44, 47, 47, 87, 89, 3, 3, 3, 567]
  nodes.each { |node| line.addNode(node) }
  nodes.uniq!
  nodes.each_index do |index| 
   assert_equal nodes[index], line.node(index), "node #{index} doesnt match" 
  end
 end

 def testLines
  lines = Lines.new
  assert_equal 0, lines.number
 end

 def testInitLines
  lines = Lines.new
  assert_equal 0, lines.number
  assert_equal EMPTY, lines.node(-1,0)
  assert_equal EMPTY, lines.node( 0,0)
  assert_equal EMPTY, lines.node( 1,0)
 end

 def testAddNodeToLines
  lines = Lines.new
  line = 1
  node = 5
  lines.addNode(line,node)
  assert_equal 2, lines.number
  assert_equal node, lines.node(line,0)
 end

 def testAddNodesToLines
  lines = Lines.new
  line = 1
  nodes = [434, 437, 223, 887, 89, 173, 4567]
  nodes.each { |node| lines.addNode(line,node) }
  nodes.each_index do |index| 
   assert_equal nodes[index], lines.node(line,index), 
       "node #{index} doesnt match" 
  end
 end

 def testRenumberLineNodes
  lines = Lines.new
  [1, 2, 3, 4].each { |node| lines.addNode(1,node) }
  [5, 6, 7].each { |node| lines.addNode(2,node) }
  lines.renumber([10,11,12,13,14,15,16,17])
  nodes1 = [11, 12, 13, 14]
  nodes2 = [15, 16, 17]
  nodes1.each_index do |index| 
   assert_equal nodes1[index], lines.node(1,index), 
       "node #{index} doesnt match" 
  end
  nodes2.each_index do |index| 
   assert_equal nodes2[index], lines.node(2,index), 
       "node #{index} doesnt match" 
  end
 end

 def testSaveLines
  lines = Lines.new
  lines.addNode(1,10)
  lines.addNode(2,20)
  lines.addNode(2,30)

  assert_equal 3, lines.number
  assert_equal EMPTY, lines.node(0,0)
  assert_equal 10, lines.node(1,0)
  assert_equal EMPTY, lines.node(1,1)
  assert_equal 20, lines.node(2,0)
  assert_equal 30, lines.node(2,1)
  assert_equal EMPTY, lines.node(2,2)

  fileName = "lineRestart"
  lines.save fileName

  restart = IO.readlines(fileName)
  assert_equal 7,  restart.size
  assert_equal 3,  restart[0].to_i
  assert_equal 0,  restart[1].to_i
  assert_equal 1,  restart[2].to_i
  assert_equal 11, restart[3].to_i
  assert_equal 2,  restart[4].to_i
  assert_equal 21, restart[5].to_i
  assert_equal 31, restart[6].to_i

  lines = Lines.new.load(fileName)
  assert_equal 3, lines.number
  assert_equal EMPTY, lines.node(0,0)
  assert_equal 10, lines.node(1,0)
  assert_equal EMPTY, lines.node(1,1)
  assert_equal 20, lines.node(2,0)
  assert_equal 30, lines.node(2,1)
  assert_equal EMPTY, lines.node(2,2)

  File.delete(fileName)
 end

end
