close all
clear all

options = optimset('jacobian','off');
MRData=rand(10);
LowerBound   = [-1 -1];
UpperBound   = [10 10];
InitialGuess = [0,1];
x = lsqnonlin(@ForwardProject,InitialGuess,LowerBound,UpperBound,options,MRData)
%MRData(1)
