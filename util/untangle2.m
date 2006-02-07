function [x,y]=untangle2(grid_x,grid_y,plot_mode)

if (nargin<2)
  grid_x = [ 0 1 1 0 ]';
  grid_y = [ 0 0 1 1 ]';
end

if (nargin<3)
  plot_mode = 0;
end

degree = size(grid_x,1);

a=[];
c=[];
for i = 1:degree
  j = i+1;
  if (i==degree)
    j = 1;
  end
  ax=-0.5*(grid_y(i)-grid_y(j));
  ay=-0.5*(grid_x(j)-grid_x(i));
  a=[a [ax;ay;1] ];
  ci = 0.5*(grid_x(i)*grid_y(j) - grid_x(j)*grid_y(i));
  c=[c;ci];
end
a = [a eye(3) ];
M=sum(c);
c = [c; [ M; M; M ]];

b=[0;0;1];

a;
c;
b;

[optx,zmin,basis] = tableau(c,a,b);

xy = a(:,basis)'\c(basis);
x = xy(1);
y = xy(2);

data = [];
for i = 1:degree
  j = i+1;
  if (i==degree)
    j = 1;
  end
  data = [ data 
	  xy(1) xy(2)
	  grid_x(i) grid_y(i) 
	  grid_x(j) grid_y(j) 
	  xy(1) xy(2)
	  ];
end

if (plot_mode)
  plot(data(:,1),data(:,2))
  disp('press key')
  pause
end

%! plot_test = 0;
%! diff_tol = 1.0e-5;

%! grid_x = [ 0 1 0 ]';
%! grid_y = [ 0 0 1 ]';
%! [x,y] = untangle2(grid_x,grid_y,plot_test);
%! assert_within(1/3,x,diff_tol)
%! assert_within(1/3,y,diff_tol)

%! grid_x = [ 0 1 1 0 ]';
%! grid_y = [ 0 0 1 1 ]';
%! [x,y] = untangle2(grid_x,grid_y,plot_test);
%! assert_within(0.5,x,diff_tol)
%! assert_within(0.5,y,diff_tol)

%! grid_x = [ 0 1.2 1.3 0 ]';
%! grid_y = [ 0 0 1.4 1.5 ]';
%! [x,y] = untangle2(grid_x,grid_y,plot_test);
%! assert_within(0.60465,x,diff_tol)
%! assert_within(0.75581,y,diff_tol)

%! grid_x = [ 0 1 0.5 0 ]';
%! grid_y = [ 0 0 0.5 1 ]';
%! [x,y] = untangle2(grid_x,grid_y,plot_test);
%! assert_within(0.25000,x,diff_tol)
%! assert_within(0.25000,y,diff_tol)

%! grid_x = [    0 1.2 0.3 0.1 ]';
%! grid_y = [ -0.1 0.0 0.2 1.5 ]';
%! [x,y] = untangle2(grid_x,grid_y,plot_test);
%! assert_within(0.11595,x,diff_tol)
%! assert_within(0.051622,y,diff_tol)

