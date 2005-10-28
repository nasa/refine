function [x,y] = data2hist(data,buckets)

if (nargin < 2 )
  buckets = 100;
end

largest  = max(data);
smallest = min(data);


delta = (largest-smallest)/buckets;

x = linspace(smallest,largest-delta,buckets)';
x = x + delta/2;

bins = data - smallest;
bins = floor(bins/delta)+1;
bins = min(bins,buckets);

y=zeros(buckets,1);

for i = 1:size(bins)
  y(bins(i)) = y(bins(i)) + 1;
end
