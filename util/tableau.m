function [optx,zmin,basis] = tableau(c,A,b)

[m,n]=size(A)

B = eye(m)
Binv = B

x = zeros(m,1)
x(m)=1

cB = c(n-m+1:n)

t = [ -cB'*b  c'-cB'*Binv*A
      Binv*b Binv*A]

basis = [n-m+1:n]


[min_reducted_cost, pivot_col] = min(t(1,2:n+1));
min_reducted_cost
pivot_col = pivot_col +1

while (min_reducted_cost<0)

  pivot_row = 0;
  step_length = inf;
  for i=2:m+1
    if (t(i,pivot_col)>0)
      this_step = t(i,1)./t(i,pivot_col)
      if (this_step<step_length)
	pivot_row=i;
	step_length = this_step;
      end
    end
  end
  pivot_row
  step_length


  disp('flip basis')
  basis(pivot_row-1) = pivot_col-1
  disp('normalize row')
  t(pivot_row,:) = t(pivot_row,:)./t(pivot_row,pivot_col)

  disp('elimate other rows')
  for i=1:m+1
    if (i!=pivot_row)
      t(i,:) = t(i,:) - t(pivot_row,:)*t(i,pivot_col);
    end
  end
  t
  [min_reducted_cost, pivot_col] = min(t(1,2:n+1));
  min_reducted_cost
  pivot_col = pivot_col +1

end

optx = zeros(n,1);
for i = 1:m
  optx(basis(i)) = t(i+1,1);
end
zmin = -t(1,1);
