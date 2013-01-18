#!/usr/bin/env ruby

raise("usage : script aflr.mapbc vgrid.mapbc") if (ARGV.size<2)

File.open(ARGV[1],'w') do |vg| 
vg.puts <<EOF
#Now
#ref.map
Patch #        BC             Family   #surf   surfIDs         Family
#---------------------------------------------------------------------
EOF
  File.open(ARGV[0]).each do |l|
    t = l.split(' ')
    next if ( t.size < 2)
    family = if (2 == t) then
               t[2]
             else
               'Addams'
             end
    patch = t[0]
    bc = case (t[1].to_i)
           when 5050
           3
           when 6662
           1
           when 3000
           5
         else
           raise "unknow bc "+t[1]
         end
    vg.printf "%d %d %d %d %s\n", patch, bc, 0, 0, family
  end
end
