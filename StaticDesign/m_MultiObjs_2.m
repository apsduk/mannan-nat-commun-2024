% ------------------------------------------------------------------------
% OBJECTIVES OF MULTIOBJECTIVE OPTIMIZATION PROBLEM 2
%
% Objectives to be minimized are:
%   - obj1 = -1 * volumetric productivity
%   - obj2 = -1 * product yield
% ------------------------------------------------------------------------

function [objs,pp] = m_MultiObjs_2(p,hPR,xPR,y0,SysTopol,maxObjVals)


% --- Redefine parameters ------------------------------------------------

% defining expression of Ep and also Tp (if expressed)...
xPR_t     = xPR;
xPR_t(12) = p(1); % ... scaling TX rate of Ep
xPR_t(13) = p(2); % ... scaling TX rate of Tp
xPR_t(10) = p(3); % ... scaling TX rate of E
xPR_t(14) = p(4); % ... scaling TX rate of TF


% --- Calculate objective values -----------------------------------------

% simulate batch culture
pp = m_SimBatchCult_ProdPerf(hPR,xPR_t,SysTopol,y0,1,0,[],[]); % pp = [vProd,pYield]

% output:
objs = (-1) * pp./maxObjVals; % this is minimised, which effectively maximises both objectives
