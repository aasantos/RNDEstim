function [hoptim,cvmat,areamat,diffmat,entropymat] = optimalbandwidth(S)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The function calculates the optimal bandwidths h_c and h_p using a
%  generalized Cross-Validation approach adapted to the situation where the
%  main aim is to estimate a risk-neutral density, associated with the
%  estimates of a second derivativa function. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: a structure S
%   S.x0          vector
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
%   S.hcmin         scalar
%   S.hcmax         scalar
%   S.hpmin         scalar
%   S.hpmax         scalar
%   S.ngrid         scalar (integer) dimension of the grid-matrix
%
%   S.sol    vector with initial estimates; it can be []
%   S.lg     string % local ; global; both
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   output: 
%      hoptim   vector   [h_c,h_p] "optimal" values
%      cvmat    matrix   sum of squared errors associated with the mean
%      areamat  matrix   mean of the areas defined for which (h_c,h_p)
%      diffmat  matrix   mean of the variation defined for which (h_c,h_p)
%      entropymat matrix mean of the entropy defined for which (h_c,h_p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
callprice = S.callprice;
callstrike = S.callstrike;
callopenint = S.callopenint;
%
putprice = S.putprice;
putstrike = S.putstrike;
putopenint = S.putopenint;
%
r = S.r;
tau = S.tau;
%
hcmin = S.hcmin;
hcmax = S.hcmax;
hpmin = S.hpmin;
hpmax = S.hpmax;
ngrid = S.ngrid;
%
hc = linspace(hcmin,hcmax,ngrid);
hp = linspace(hpmin,hpmax,ngrid);
%
cvmat = zeros(ngrid,ngrid);
areamat = zeros(ngrid,ngrid);
diffmat = zeros(ngrid,ngrid);
entropymat = zeros(ngrid,ngrid);
%
strikes = [callstrike;putstrike];
strikeunique = unique(strikes);
%
%%% defines the range of the strikes to be considered. It uses the quantile
%%% function to exclude the extreme strike values, but other form of
%%% excluding extreme strike values can be used
strikeunique = strikeunique( (strikeunique > quantile(strikes,0.01)) &  ...
                             (strikeunique < quantile(strikes,0.99))); 
nstrike = length(strikeunique);
%
%
%%% defines the number of points used to approximate the rnd for each
%%% iterations and through it calculate area, variation and entropy. This
%%% value can be changed to define a grid with less or more grid points
m = 100;
%
%
jj = 0;
for k=1:ngrid
    for j=1:ngrid
        jj = jj + 1;
        fprintf("%d/%d\n",jj,ngrid*ngrid);
        cvv = 0.0;
        area = 0.0;
        dd = 0.0;
        ent = 0.0;
        for i=1:nstrike
            %%% defines the range using the quantile function to define the
            %%% lower and upper limits for the m points. The range can be
            %%% changed by changing the quantiles which define the limits
            x0 = linspace(quantile(strikes,0.025),quantile(strikes,0.975),m);
            xx0 = strikeunique(i);
            x0 = sort([x0 xx0]);
            index = find(x0 == xx0);
            index = index(1);
            %
            callstriket = callstrike(callstrike ~= xx0);
            callpricet = callprice(callstrike ~= xx0);
            callopenintt =  callopenint(callstrike ~= xx0);
            putstriket = putstrike(putstrike ~= xx0);
            putpricet = putprice(putstrike ~= xx0);
            putopenintt = putopenint(putstrike ~= xx0);
            %
            [fcall,fput,~,ddf] = callputfunctionf(callstriket,callpricet, ...
                callopenintt,putstriket,putpricet,putopenintt,x0,hc(k),hp(j),r,tau);
            %
            area = area + areadensity(x0,ddf);
            dd = dd + sum(abs(diff(ddf,1)));
            ent = ent + entropy(x0,ddf);
            call = callprice(callstrike==xx0);
            put = putprice(putstrike==xx0);
            nc = length(call);
            np = length(put);
            cvc = 0.0;
            cvp = 0.0;
            if nc > 0
                cvc = (1/nc)*sum((call - fcall(index)).^2);
            end
            if np > 0
                cvp = (1/np)*sum((put - fput(index)).^2);
            end
            cvv = cvv + (cvc + cvp);
        end
        cvv = cvv/nstrike;
        area = area/nstrike;
        dd = dd/nstrike;
        ent = ent/nstrike;
        cvmat(k,j) = cvv;
        areamat(k,j) = area;
        diffmat(k,j) = dd;
        entropymat(k,j) = ent;
    end
end
matcrit = cvmat.*diffmat + (1 + abs(areamat -1))./entropymat;
[row,col] = minMatrix(matcrit);
hcop = hc(row);
hpop = hp(col);
hoptim = [hcop;hpop];
%
end