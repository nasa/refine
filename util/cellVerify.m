

# cat fort.2?? | sort > allSorted
# sed -e 's/[a-z]/ /g' < allSorted > clean

if (exist('clean')!=1)
 load clean
endif

cell_col = 1
local_node_1 = 8
global_node1 = 4
nnodes0 = 12;
length = size(clean,1)

local = (clean(:,local_node_1)<=clean(:,nnodes0));
local_test = sum(local)
local_cells = find(local);

max_cell = max(clean(:,cell_col))
min_cell = min(clean(:,cell_col))

hits = zeros(max_cell,1);
hits(clean(local_cells,cell_col))++;

[max_hits, max_hits_index] = max(hits)
[min_hits, min_hits_index] = min(hits)
hits--;

bad_cells = find(hits)

hits(bad_cells)

nodes = zeros(4,max_cell);
for line = 1:max_cell
  cell = clean(line,cell_col);
  if (nodes(1,cell)==0)
    nodes(1,cell) = clean(line,global_node1+0);
    nodes(2,cell) = clean(line,global_node1+1);
    nodes(3,cell) = clean(line,global_node1+2);
    nodes(4,cell) = clean(line,global_node1+3);
  else
    if (nodes(1,cell) != clean(line,global_node1+0) || \
	nodes(2,cell) != clean(line,global_node1+1) || \
	nodes(3,cell) != clean(line,global_node1+2) || \
	nodes(4,cell) != clean(line,global_node1+3) ) 
      line
      cell
      nodes(:,cell)
      clean(line,global_node1:global_node1+3)
      clean(line,:)
    endif
  endif
endfor

				# line 283532
				# cell 251211
				# first 3952 4351 4222 136
				# next 32065  57518  46704  18643
