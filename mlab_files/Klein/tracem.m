function y = tracem(x);

n = size(x,2);
m = size(x,1)/n;
for i=1:m
   y(i,1)=trace(x((n*(i-1)+1):i*n,1:n));
end
