path(path,'../util');

histogram;

[a,b] = data2hist(edge_length);
[x,y] = bar(a,b);
xy=[x y];

gset term postscript eps enhanced
gset output 'histogram_edge_length.eps'
gset xlabel 'Segment Length in Metric'
gset ylabel 'Incidence'
gplot [:] [:] xy with lines 01 title ''
closeplot

[a,b] = data2hist(triangle_uv_area);
[x,y] = bar(a,b);
xy=[x y];

gset term postscript eps enhanced
gset output 'histogram_triangle_uv_area.eps'
gset xlabel 'Trangle UV Area'
gset ylabel 'Incidence'
gplot [:] [:] xy with lines 01 title ''
closeplot

[a,b] = data2hist(tet_volume);
[x,y] = bar(a,b);
xy=[x y];

gset term postscript eps enhanced
gset output 'histogram_tet_volume.eps'
gset xlabel 'Tetrahedral Volume'
gset ylabel 'Incidence'
gplot [:] [:] xy with lines 01 title ''
closeplot
