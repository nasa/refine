

# cat fort.2?? | sort > allSorted
# sed -e 's/[a-z]/ /g' < allSorted > clean

if (exist('clean')!=1)
 load clean
endif

cell_col = 1
local_node_1 = 8
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

