function SO = npcallputoptimLG(S)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The function estimate the call pricing function, the put pricing function
%  and the risk-neutral density. It can use two methods, a sequential
%  approach that evaluate the value of the function at each point on a
%  previous defined grid for x (strikes). The gobal approach estimates the
%  functions in one go. The first approach, assuming that the vector x has
%  n elements, solves n Quadratic Programming problems with a decision
%  vector of dimension 8. The second approach compacts all in one problem
%  with a decision vector of dimension n*8, allowing additional constraints
%  to be imposed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: 
%   A structure S
%   S.x0          vector  (support defined by a grid of points)
%   S.putprice    vector
%   S.putstrike   vector
%   S.putopenint  vector
%
%   S.callprice   vector
%   S.callstrike  vector
%   S.callopenint vector
%
%   S.r         scalar (annual interest rate)
%   S.tau       scalar (time to maturity)
%
%   S.hc        scalar (bandwidth for calls)
%   S.hp        scalar (bandwidth for puts)
%
%   S.sol    vector  (initial estimates; it can be [])
%   S.lg     string  (type of optimiation: "local" ; "global" ; "both")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   output:
%   A structure SO
%    
%    SO.x     (supported used)
%    SO.call  (call mean pricing function)
%    SO.put   (put mean pricing function)
%    SO.dcall (call first defivative function; put obtained by adding +1)
%    SO.ddcall (second derivative function; coincide for calls and puts)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function z0 = kernormal(x)
        z0 = (1/sqrt(2*pi))*exp(-0.5*x.^2);
    end
%
%
    function [A,b,Aeq,beq,lb,ub] = restlocal(x0)
        %
        m1 = length(x0);
        %
        A = [];
        b = [];
        Aeq = zeros(m1*3,8*m1);
        %         % impose that the third derivatives are equal
        col = 4;
        for i=1:m1
            Aeq(i,col) = 0.0;
            Aeq(i,col+4) = 0.0;
            col = col + 8;
        end
        % impose that the second derivatives are equal
        col = 3;
        for i=1:m1
            Aeq(i+m1,col) = 1.0;
            Aeq(i+m1,col+4) = -1.0;
            col = col + 8;
        end
        % impose that b2 and b6 (first derivatives) b2 = b6 - 1
        col = 2;
        for i=1:m1
            Aeq(i+2*m1,col) = -1.0;
            Aeq(i+2*m1,col + 4) = 1;
            col = col + 8;
        end
        beq = [zeros(m1,1);zeros(m1,1);exp(-r*tau)*ones(m1,1)];
        %
        Aeq = sparse(Aeq);
        %
        lbt = [0       -exp(-r*tau)      0   -Inf       0    0               0   -Inf];
        ubt = [Inf          0          Inf    Inf      Inf   exp(-r*tau)    Inf   Inf];
        lb = repmat(lbt,1,m1);
        ub = repmat(ubt,1,m1);
    end
%
%
%
    function [A,b,Aeq,beq,lb,ub] = restglobal(x0)
        %%%
        function rr = sumw(xx0)
            nn = length(xx0);
            rr = zeros(1,nn);
            for ii=1:nn
                if ii==1
                    rr(ii) = 0.5*(xx0(ii+1) - xx0(ii));
                elseif ii==nn
                    rr(ii) = 0.5*(xx0(ii) - xx0(ii-1));
                else
                    rr(ii) = 0.5*(xx0(ii+1) - xx0(ii-1));
                end
            end
        end
        %%%
        m2 = length(x0);
        k = m2 - 1;
        %%%
        Aeq = zeros(m2*3,8*m2);
        col = 4;
        for i=1:m2
            Aeq(i,col) = 0.0;
            Aeq(i,col+4) = 0.0;
            col = col + 8;
        end
        col = 3;
        for i=1:m2
            Aeq(i+m2,col) = 1.0;
            Aeq(i+m2,col+4) = -1.0;
            col = col + 8;
        end
        col = 2;
        for i=1:m2
            Aeq(i+2*m2,col) = -1.0;
            Aeq(i+2*m2,col + 4) = 1;
            col = col + 8;
        end
        beq = [zeros(m2,1);zeros(m2,1);exp(-r*tau)*ones(m2,1)];
        %%%
        vv = sumw(x0);
        Aeq1 = zeros(1,8*m2);
        Aeq1(3:8:(8*m2)) = vv;
        Aeq = [Aeq;Aeq1];
        beq = [beq;exp(-r*tau)];
        %%%
        A = zeros(2*k,8*m2);
        col  = 1;
        for i=1:k
            A(i,col) = -1.0;
            A(i,col+8) = 1.0;
            col = col + 8;
        end
        col = 5;
        for i=1:k
            A(i+k,col) = 1.0;
            A(i+k,col+8) = -1.0;
            col = col + 8;
        end
        b = zeros(2*k,1);
        %%%
        A3 = zeros(2*k,8*m2);
        col = 2;
        for i=1:k
            A3(i,col) = 1.0;
            A3(i,col+8) = -1.0;
            col = col + 8;
        end
        col = 6;
        for i=1:k
            A3(i+k,col) = 1.0;
            A3(i+k,col+8) = -1.0;
            col = col + 8;
        end
        A = [A;A3];
        b = [b;zeros(2*k,1)];
        %%%
        A = sparse(A);
        Aeq = sparse(Aeq);
        %%%
        lbt = [0       -exp(-r*tau)     0    -Inf       0        0         0    -Inf];
        ubt = [Inf          0          Inf    Inf       Inf   exp(-r*tau)  Inf   Inf];
        lb = repmat(lbt,1,m2);
        ub = repmat(ubt,1,m2);
    end
%
%
    function [H,f] = buildHf(x0)
        m3 = length(x0);
        f = zeros(8*m3,1);
        row = [];
        col = [];
        for i=1:2*m3
            l1 = (i-1)*4 + 1;
            l2 = l1 + 3;
            [Xcol,Yrow] = meshgrid(l1:l2,l1:l2);
            row = [row;Yrow(:)];
            col = [col;Xcol(:)];
        end
        elem = zeros(32*m3,1);
        for i=1:m3
            %
            w1 = sqrt(callopenint.*kernormal((callstrike - x0(i))/hc));
            yc = w1.*callprice;
            x1 = w1;
            x2 = x1.*(callstrike - x0(i));
            x3 = (1/2)*x2.*(callstrike - x0(i));
            x4 = (2/6)*x3.*(callstrike - x0(i));
            X = [x1 x2 x3 x4];
            XXc = X'*X;
            Xyc = X'*yc;
            w1 = sqrt(putopenint.*kernormal((putstrike - x0(i))/hp));
            yp = w1.*putprice;
            x1 = w1;
            x2 = x1.*(putstrike - x0(i));
            x3 = (1/2)*x2.*(putstrike - x0(i));
            x4 = (2/6)*x3.*(putstrike - x0(i));
            X = [x1 x2 x3 x4];
            XXp = X'*X;
            Xyp = X'*yp;
            y = [Xyc;Xyp];
            index = ((i-1)*8+1):(i*8);
            f(index) = y;
            index1 = ((i-1)*32 +1):(i*32);
            elem(index1) = [XXc(:);XXp(:)];
        end
        H = sparse(row,col,elem,8*m3,8*m3);
        f = -1.0*f;
    end
%
x0 = S.x0;
%
callstrike = S.callstrike;
callprice = S.callprice;
callopenint = S.callopenint;
%
putstrike = S.putstrike;
putprice = S.putprice;
putopenint = S.putopenint;
%
r = S.r;
tau = S.tau;
%
hp = S.hp;
hc = S.hc;
sol = S.sol;
%
m = length(x0);
[H,f] = buildHf(x0);
%
%
options = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','off');
%
%
if strcmp("global",S.lg)
        [A,b,Aeq,beq,lb,ub] = restglobal(x0);
    out = quadprog(H,f,A,b,Aeq,beq,lb,ub,sol,options);
    call = zeros(m,1);
    put = zeros(m,1);
    dcall = zeros(m,1);
    ddcall = zeros(m,1);
    for iter=1:m
        ind1 = (iter-1)*8 + 1;
        ind5 = (iter-1)*8 + 5;
        ind2 = (iter-1)*8 + 2;
        ind3 = (iter-1)*8 + 3;
        call(iter) = out(ind1);
        put(iter) = out(ind5);
        dcall(iter) = out(ind2);
        ddcall(iter) = out(ind3);
    end
    %
    SO.callG = call;
    SO.putG = put;
    SO.dcallG = dcall;
    SO.ddcallG = ddcall;
    SO.sol = out;
    %
    %
    [A,b,Aeq,beq,lb,ub] = restglobal(x0);
    out = quadprog(H,f,A,b,Aeq,beq,lb,ub,sol,options);
    call = zeros(m,1);
    put = zeros(m,1);
    dcall = zeros(m,1);
    ddcall = zeros(m,1);
    for iter=1:m
        ind1 = (iter-1)*8 + 1;
        ind5 = (iter-1)*8 + 5;
        ind2 = (iter-1)*8 + 2;
        ind3 = (iter-1)*8 + 3;
        call(iter) = out(ind1);
        put(iter) = out(ind5);
        dcall(iter) = out(ind2);
        ddcall(iter) = out(ind3);
    end
    %
    SO.callG = call;
    SO.putG = put;
    SO.dcallG = dcall;
    SO.ddcallG = ddcall;
    SO.sol = out;
end


if strcmp("local",S.lg)
    [A,b,Aeq,beq,lb,ub] = restlocal(x0);
    out = quadprog(H,f,A,b,Aeq,beq,lb,ub,sol,options);
    call = zeros(m,1);
    put = zeros(m,1);
    dcall = zeros(m,1);
    ddcall = zeros(m,1);
    for iter=1:m
        ind1 = (iter-1)*8 + 1;
        ind5 = (iter-1)*8 + 5;
        ind2 = (iter-1)*8 + 2;
        ind3 = (iter-1)*8 + 3;
        call(iter) = out(ind1);
        put(iter) = out(ind5);
        dcall(iter) = out(ind2);
        ddcall(iter) = out(ind3);
    end
    %
    SO.x = x0;
    SO.call = call;
    SO.put = put;
    SO.dcall = dcall;
    SO.ddcall = ddcall;     
    SO.sol = out;
end
%
%
%
if strcmp("both",S.lg)
    [A,b,Aeq,beq,lb,ub] = restglobal(x0);
    out = quadprog(H,f,A,b,Aeq,beq,lb,ub,sol,options);
    call = zeros(m,1);
    put = zeros(m,1);
    dcall = zeros(m,1);
    ddcall = zeros(m,1);
    for iter=1:m
        ind1 = (iter-1)*8 + 1;
        ind5 = (iter-1)*8 + 5;
        ind2 = (iter-1)*8 + 2;
        ind3 = (iter-1)*8 + 3;
        call(iter) = out(ind1);
        put(iter) = out(ind5);
        dcall(iter) = out(ind2);
        ddcall(iter) = out(ind3);
    end
    %
    SO.callG = call;
    SO.putG = put;
    SO.dcallG = dcall;
    SO.ddcallG = ddcall;
    SO.sol = out;
    %
    %
    [A,b,Aeq,beq,lb,ub] = restglobal(x0);
    out = quadprog(H,f,A,b,Aeq,beq,lb,ub,sol,options);
    call = zeros(m,1);
    put = zeros(m,1);
    dcall = zeros(m,1);
    ddcall = zeros(m,1);
    for iter=1:m
        ind1 = (iter-1)*8 + 1;
        ind5 = (iter-1)*8 + 5;
        ind2 = (iter-1)*8 + 2;
        ind3 = (iter-1)*8 + 3;
        call(iter) = out(ind1);
        put(iter) = out(ind5);
        dcall(iter) = out(ind2);
        ddcall(iter) = out(ind3);
    end
    %
    SO.callG = call;
    SO.putG = put;
    SO.dcallG = dcall;
    SO.ddcallG = ddcall;
    SO.sol = out;
    %
    %
    [A,b,Aeq,beq,lb,ub] = restlocal(x0);
    out = quadprog(H,f,A,b,Aeq,beq,lb,ub,sol,options);
    call = zeros(m,1);
    put = zeros(m,1);
    dcall = zeros(m,1);
    ddcall = zeros(m,1);
    for iter=1:m
        ind1 = (iter-1)*8 + 1;
        ind5 = (iter-1)*8 + 5;
        ind2 = (iter-1)*8 + 2;
        ind3 = (iter-1)*8 + 3;
        call(iter) = out(ind1);
        put(iter) = out(ind5);
        dcall(iter) = out(ind2);
        ddcall(iter) = out(ind3);
    end
    %
    SO.x = x0;
    SO.call = call;
    SO.put = put;
    SO.dcall = dcall;
    SO.ddcall = ddcall;    
    %
end
%
%
end