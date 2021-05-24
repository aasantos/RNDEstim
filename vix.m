%%
clc;clear;
load vixData.mat callstrike callprice callopenint ...
                 putstrike putprice putopenint
%%
%
index = callstrike > 9 & callstrike < 41;
callstrike = callstrike(index);
callprice = callprice(index);
callopenint = callopenint(index);
index = putstrike > 9 & putstrike < 41;
putstrike = putstrike(index);
putprice = putprice(index);
putopenint = putopenint(index);
%
%
%% Optimal bandwith choice (It may take several hours to run)
%
%  It uses the function optimalbandwidth function, which has as input a 
%  structure. The function outputs a vector with the optimal values for the
%  h_c and h_p, and the respective matrices used to calculate such values.
%
S.callprice = callprice;
S.callstrike = callstrike;
S.callopenint = callopenint;
%
S.putprice = putprice;
S.putstrike = putstrike;
S.putopenint = putopenint;
%
S.r = 0.02 ;
S.tau = 1/12;
%
S.hcmin = 1.0;
S.hcmax = 2.0;
S.hpmin = 1.0;
S.hpmax = 2.0;
S.ngrid = 53;
%
[hoptim,cvmat,areamat,diffmat,entropymat] = optimalbandwidth(S);
%% Bootstrap Analysis 
%
%  It uses the function cibootstrap, which has as input a structure, where
%  named SI with the corresponding fields. The function gives a structured
%  where are saved as fields, for each iteration, using a sequential local
%  and a global approach the function's mean values, first derivatives and
%  second derivatives. The distribution of such objects can then be
%  approximated by simulation obtained through bootstrap
%
x0 = linspace(10.2,37,300);
r = 0.02 ;
tau = 1/12;
%
niter = 5000;
SI.callstrike = callstrike;
SI.callprice = callprice;
SI.callopenint = callopenint;
SI.putstrike = putstrike;
SI.putprice = putprice;
SI.putopenint = putopenint;
SI.r = r;
SI.tau = tau;
SI.x0 = x0;
SI.hcop = hoptim(1);
SI.hpop = hoptim(2);
SI.niter = niter;
%
SB = cibootstrap(SI);
%% Graphical analysis
%
%  comparison of sample and estimated mean function for calls and puts
%  using the sequential local and the global approach
%
x0 = linspace(10.2,37,100);
S.callstrike = callstrike;
S.callprice = callprice;
S.callopenint = callopenint;
S.putstrike = putstrike;
S.putprice = putprice;
S.putopenint = putopenint;
S.r = r;
S.tau = tau;
S.sol = [];
S.lg = "both";
S.x0 = x0;
S.hc = hoptim(1);
S.hp = hoptim(2);
%
S1 = npcallputoptimLG(S);
%
subplot(1,2,1)
plot(callstrike,callprice,'.','color','blue')
xlabel('strikes','FontSize',14)
ylabel('prices','FontSize',14)
title('S&P500 - call and put; sample and mean','FontSize',18)
hold on
plot(x0,S1.call,'+','color','blue')
plot(S1.x,S1.callG,'o','color','red')
xlim([23.5 28.25])
plot(putstrike,putprice,'.','color','red')
plot(x0,S1.put,'+','color','blue')
plot(S1.x,S1.putG,'o','color','red')
hold off
%  
%  Comparison of the risk-neutral densities
%
x0 = linspace(10.2,37,300);
S.callstrike = callstrike;
S.callprice = callprice;
S.callopenint = callopenint;
S.putstrike = putstrike;
S.putprice = putprice;
S.putopenint = putopenint;
S.r = r;
S.tau = tau;
S.x0 = x0;
S.hc = hoptim(1);
S.hp = hoptim(2);
S.sol = [];
S.lg = "both";
%
S1 = npcallputoptimLG(S);
%
subplot(1,2,2)
plot(x0,exp(r*tau)*S1.ddcall,'LineWidth',1.25,'color','blue')
xlabel('S_T','FontSize',14)
ylabel('density','FontSize',14)
title('S&P500 - risk-neutral density','FontSize',18)
hold on
plot(x0,exp(r*tau)*S1.ddcallG,'color','red','LineWidth',1.25)
xlim([10.2 37])
hold off
%% Results from the bootstrap analysis
%
%
x0 = linspace(10.2,37,300);
rndsample = SB.rndsample;
rndsampleG = SB.rndsampleG;
nxy = length(x0);
qqval = [];
qqvalg = [];
for i=1:nxy
    qq = quantile(rndsample(i,:),[0.05 0.5 0.95]);
    qqval = [qqval;qq];
    qq = quantile(rndsampleG(i,:),[0.05 0.5 0.95]);
    qqvalg = [qqvalg;qq];
end
%
subplot(1,2,1)
plot(x0,qqval(:,3),'--','color','blue','LineWidth',1.25)
xlabel('S_T','FontSize',12)
ylabel('density','FontSize',12)
title('S&P500 - risk-neutral density quantiles','FontSize',16)
xlim([10.2 37])
hold on
plot(x0,qqvalg(:,3),'--','color','red','LineWidth',1.25)
%
plot(x0,qqval(:,2),'color','blue','LineWidth',1.25)
plot(x0,qqvalg(:,2),'color','red','LineWidth',1.25)
%
plot(x0,qqval(:,1),'--','color','blue','LineWidth',1.25)
plot(x0,qqvalg(:,1),'--','color','red','LineWidth',1.25)
hold off
%
%
%  Difference of results in the bootstrap analysis
%
diff1 = qqval(:,3) - qqval(:,1);
diff2 = qqvalg(:,3) - qqvalg(:,1);
%
subplot(1,2,2)
plot(x0,diff1,'LineWidth',1.25,'color','blue')
xlabel('S_T','FontSize',14)
ylabel('difference','FontSize',14)
title('S&P500 - difference of quantiles','FontSize',18)
xlim([10.2 37])
hold on
plot(x0,diff2,'LineWidth',1.25,'color','red')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
