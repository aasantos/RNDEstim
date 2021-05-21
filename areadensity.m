function out = areadensity(x0,y)
        n = length(x0);
        rr = zeros(1,n);
        for i=1:n
            if i==1
                rr(i) = 0.5*(x0(i+1) - x0(i))*y(i);
            elseif i==n
                rr(i) = 0.5*(x0(i) - x0(i-1))*y(i);
            else
                rr(i) = 0.5*(x0(i+1) - x0(i-1))*y(i);
            end
        end
        out = sum(rr);
end