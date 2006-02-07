
grid_x = [ 0 1 1 0 ]';
grid_y = [ 0 0 1 1 ]';
#grid_x = [ 0 1.2 1.3 0 ]';
#grid_y = [ 0 0 1.4 1.5 ]';
#grid_x = [ 0 1 0.5 0 ]';
#grid_y = [ 0 0 0.5 1 ]';
#grid_x = [    0 1.2 0.3 0.1 ]';
#grid_y = [ -0.1 0.0 0.2 1.5 ]';

degree = size(grid_x,1)

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

a
c
b

[optx,zmin,basis] = tableau(c,a,b)

xy = a(:,basis)'\c(basis)


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

if (1)
  plot(data(:,1),data(:,2))
  disp('press key')
  pause
end

