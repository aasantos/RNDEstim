function SB = cibootstrap(SI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This function implements a bootstrap analysis to define the confidence
% intervals associated with the Risk-Neutral density estimation using
% nonparametric methods.  
%
%  input:
%  This function as one input in the form of a
%  structure SI, and the fields of the structure must be:
%     - callstrike:
%     - callprice:
%     - callopenint: this variable means weights, openinterest or volume
%
%     - putstrike:
%     - putprice:
%     - putopenint: this variable means weights, openinterest or volume
%
%     - x0:  the range (grid points used to evaluate the functions)
%     - r:   risk-free interest rate
%     - tau: time to maturity in years
%
%     - hcop: bandwidth for calls
%     - hpop: bandwidth for puts
%
%     - niter: number of iterations used in the bootstrap simulation
%
%
%     ouput:
%     The function outputs a structure with the mean functions for calls
%     and puts, the first derivative and the second derivative. For these
%     three functions, the sequential and the global (G) version is
%     recorded. The fields for the SB structure are
%   
%      - callmeansample
%      - callmeansampleG
%      - putmeansample
%      - putmeansampleG
%      
%      - callfirstderivative
%      - callfirstderivativeG: the put first derivative differs only by a
%      constant)
%      - rndsample
%      - rndsampleG
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
callstrike = SI.callstrike;
callprice = SI.callprice;
callopenint = SI.callopenint;
%
putstrike = SI.putstrike;
putprice = SI.putprice;
putopenint = SI.putopenint;
%
x0 = SI.x0;
r = SI.r;
tau = SI.tau;
%
hcop = SI.hcop;
hpop = SI.hpop;
%
niter = SI.niter;
%
ncalls = length(callprice);
nputs = length(putprice);
callindex = 1:ncalls;
putindex = 1:nputs;
weigcall = callopenint/sum(callopenint);
weigput = putopenint/sum(putopenint);
nxy = length(x0);
%
%
callmeansample = zeros(nxy,niter);
putmeansample = zeros(nxy,niter);
callfirstderiv = zeros(nxy,niter);
rndsample = zeros(nxy,niter);
%
callmeansampleG = zeros(nxy,niter);
putmeansampleG = zeros(nxy,niter);
callfirstderivG = zeros(nxy,niter);
rndsampleG = zeros(nxy,niter);
%
S.r = r;
S.tau = tau;
S.hc = hcop;
S.hp = hpop;
S.x0 = x0;
S.sol = [];
S.lg = "both";
%
h = waitbar(0,'Please wait ... ');
for i=1:niter
    waitbar(i/niter,h,'Please wait ....');
    %
    indexc = datasample(callindex,ncalls,'Replace',true,'Weights',weigcall);
    indexp = datasample(putindex,nputs,'Replace',true,'Weights',weigput);
    %
    callpricet = callprice(indexc);
    callstriket = callstrike(indexc);
    callopenintt = callopenint(indexc);
    %
    putpricet = putprice(indexp);
    putstriket = putstrike(indexp);
    putopenintt = putopenint(indexp);
    %
    S.callprice = callpricet;
    S.callstrike = callstriket;
    S.callopenint = callopenintt;
    S.putprice = putpricet;
    S.putstrike = putstriket;
    S.putopenint = putopenintt;
    %
    SO = npcallputoptimLG(S);
    S.sol = SO.sol;
    %
    callmeansample(:,i) = SO.call;
    putmeansample(:,i) = SO.put;
    callfirstderiv(:,i) = SO.dcall;
    rndsample(:,i) = exp(r*tau)*SO.ddcall;
    %
    callmeansampleG(:,i) = SO.callG;
    putmeansampleG(:,i) = SO.putG;
    callfirstderivG(:,i) = SO.dcallG;
    rndsampleG(:,i) = exp(r*tau)*SO.ddcallG;
    %
    % This section can be uncommented to debug purposes, and after
    % iteration 50 it shows the evolution of the confidence intervals at
    % each new iteration performed
    %
%     if i > 50
%         nxy = length(x0);
%         qqval = [];
%         qqvalg = [];
%         for j=1:nxy
%             qq = quantile(rndsample(j,1:i),[0.1 0.5 0.9]);
%             qqval = [qqval;qq];
%             qq = quantile(rndsampleG(j,1:i),[0.1 0.5 0.9]);
%             qqvalg = [qqvalg;qq];
%         end
%         %
%         plot(x0,qqval(:,3),'--','color','blue')
%         hold on
%         plot(x0,qqvalg(:,3),'--','color','red')
%         %
%         plot(x0,qqval(:,2),'color','blue')
%         plot(x0,qqvalg(:,2),'color','red')
%         %
%         plot(x0,qqval(:,1),'--','color','blue')
%         plot(x0,qqvalg(:,1),'--','color','red')
%         drawnow
%         hold off
%     end
    %
    %
end
close(h);
%
%
SB.callmeansample = callmeansample;
SB.callmeansampleG = callmeansampleG;
SB.putmeansample = putmeansample;
SB.putmeansampleG = putmeansampleG;
SB.callfirstderiv = callfirstderiv;
SB.callfirstderivG = callfirstderivG;
SB.rndsample = rndsample;
SB.rndsampleG = rndsampleG;
%
%
end