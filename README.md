# RNDEstim
Estimate RND functions using option prices
Data available: sample of data from VIX and S&P500 options contracts, from April 16 to April 20, 2018, 
with maturity date at May 18, 2018. The MATLAB data format files are vixData.mat and sp500Data.mat
There is a series of functions used to implement the RND estimation and respective analysis:

areadensity.m

entropy.m

minMatrix.m

callputfunctionf.m

npcallputoptimLG.m

optimalbandwidth.m

cibootstrap.m

#

Two scripts (functions) are defined applied to do the analysis for VIX and S&P500
vix.m
sp500.m

#
By running the scripts (functions), it loads the data, calculates the optimal
bandwidth values, estimates the respective RNDs, sequential local and global, and
implements a bootstrap analysis to define confidence intervals. The results are
compared with a series of graphic analysis.
