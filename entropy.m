function out = entropy(x0,y0)
n = length(x0);
rr = zeros(1,n);
a = areadensity(x0,y0);
y = y0./a;
for i=1:n
    if i==1
            rr(i) = 0.5*(x0(i+1) - x0(i))*y(i)*log(y(i) + eps);
    elseif i==n
            rr(i) = 0.5*(x0(i) - x0(i-1))*y(i)*log(y(i) + eps);
    else
            rr(i) = 0.5*(x0(i+1) - x0(i-1))*y(i)*log(y(i) + eps);
    end
end
out = -1.0*sum(rr);
end