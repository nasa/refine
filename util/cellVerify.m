

# cat fort.2?? | sort | sed -e 's/[a-z]/ /g' > clean

if (exist('clean')!=1)
 load clean
endif

cell_col =1
length = size(clean,1)

if (1==0)
  cell = clean(:,cell_col);
  shift = [clean(2:length,cell_col);clean(1,cell_col)];
  dup = (cell==shift);

  for i=1:length;
    if (dup(i))
      printf( "%6d%6d%6d%6d%6d %6d%6d p%3d\n",
	     clean(i,cell_col),clean(i,4:7),clean(i,8),clean(i,12),clean(i,2));
    endif
    if (dup(i))
      printf( "%6d%6d%6d%6d%6d %6d%6d p%3d\n",
	     clean(i+1,cell_col),clean(i+1,4:7),clean(i+1,8),clean(i+1,12),clean(i+1,2));
    endif
  endfor
endif

local = (clean(:,8)<=clean(:,12));
local_test = sum(local);
local_cells = find(local)

max_cell = max(clean(:,cell_col))
min_cell = min(clean(:,cell_col))


hits = zeros(max_cell,1);
hits(clean(local_cells,cell_col))++;
hits--;

bad_cells = find(hits)
