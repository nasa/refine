function [optx,zmin,basis] = tableau(c,A,b)

[m,n]=size(A);

B = eye(m);
Binv = B;

x = zeros(m,1);
x(m)=1;

cB = c(n-m+1:n);

t = [ -cB'*b  c'-cB'*Binv*A
      Binv*b Binv*A];

basis = [n-m+1:n];
in_basis = zeros(1,n);
in_basis(basis) = 1:m;

min_reducted_cost = -1.0;
while (min_reducted_cost<0)

  min_reducted_cost = 1.0;
  best_divisor = 0.0;
  for j = 2:n+1
    if (0==in_basis(j-1))
      if ( t(1,j) < 0 )
	pivot = 0;
	step_length = inf;
	for i=2:m+1
	  if (t(i,j)>0)
	    this_step = t(i,1)./t(i,j);
	    if (this_step<step_length)
	      pivot=i;
	      step_length = this_step;
	      divisor = t(i,j);
	    end
	  end
	end
	if ( 0 != pivot && abs(divisor) > abs(best_divisor))
	  pivot_row = pivot;
	  pivot_col = j;
	  min_reducted_cost = t(1,pivot_col);
	  best_divisor = divisor;
	end
      end
    end
  end

  if (min_reducted_cost < 0)
    best_divisor;
    in_basis(basis(pivot_row-1)) = 0;
    basis(pivot_row-1) = pivot_col-1;
    in_basis(pivot_col-1) = pivot_row-1;

    t(pivot_row,:) = t(pivot_row,:)./t(pivot_row,pivot_col);
    for i=1:m+1
      if (i!=pivot_row)
	t(i,:) = t(i,:) - t(pivot_row,:)*t(i,pivot_col);
      end
    end
    t;
  end
end

optx = zeros(n,1);
for i = 1:m
  optx(basis(i)) = t(i+1,1);
end
zmin = -t(1,1);
