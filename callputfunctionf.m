function [fcall,fput,df,ddf] = callputfunctionf(callstrike,callprice,callopenint, ...
    putstrike,putprice,putopenint,x0,hc,hp,r,tau)
    %
    %
    function [y,X] = XX(strike,price,openint,x0,h)
        %
        function z0 = kernormal(x)
            z0 = (1/sqrt(2*pi))*exp(-0.5*x.^2);
        end
        %
        w = sqrt(openint.*kernormal((strike - x0)/h));
        y = w.*price;
        x1 = w;
        x2 = x1.*(strike - x0);
        x3 = (1/2)*x2.*(strike - x0);
        x4 = (2/6)*x3.*(strike - x0);
        X = [x1 x2 x3 x4];
        %
        %
    end
    %
    %
m = length(x0);
fcall = zeros(m,1);
fput = zeros(m,1);
df = zeros(m,1);
ddf = zeros(m,1);
%
nc = length(callstrike);
np = length(putstrike);
options = optimoptions('lsqlin','Display','off');
%
Aeq = [0 0  -1  0 0 0 1 0;
    0 -1  0  0 0 1 0 0];
beq = [0;exp(-r*tau)];
lb = [0   -exp(-r*tau) 0  -Inf   0       0         0   -Inf];
ub = [Inf  0           Inf Inf   Inf exp(-r*tau)  Inf  Inf];

[yc,XXc] = XX(callstrike,callprice,callopenint,x0(1),hc);
[yp,XXp] = XX(putstrike,putprice,putopenint,x0(1),hp);

X = [XXc zeros(nc,4);zeros(np,4) XXp];
y = [yc;yp];
beta = lsqlin(X,y,[],[],Aeq,beq,lb,ub,[],options);
fcall(1) = beta(1);
fput(1) = beta(5);
df(1) = beta(2);
ddf(1) = beta(3);

for i=2:m
    [yc,XXc] = XX(callstrike,callprice,callopenint,x0(i),hc);
    [yp,XXp] = XX(putstrike,putprice,putopenint,x0(i),hp);
    X = [XXc zeros(nc,4);zeros(np,4) XXp];
    y = [yc;yp];
    beta = lsqlin(X,y,[],[],Aeq,beq,lb,ub,beta,options);
    fcall(i) = beta(1);
    fput(i) = beta(5);
    df(i) = beta(2);
    ddf(i) = beta(3);
end
%
end
