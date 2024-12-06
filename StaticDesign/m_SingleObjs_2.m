% ------------------------------------------------------------------------
% COMPUTING SELECTED SINGLE OBJECTIVE OF MULTIOBJECTIVE OPTIM. PROBLEM 2
%
% Objective to be minimized is either:
%   - obj1 = -1 * volumetric productivity
%   - obj2 = -1 * product yield
% ------------------------------------------------------------------------

function obj = m_SingleObjs_2(p,hPR,xPR,y0,SysTopol,objID)


% --- Redefine parameters ------------------------------------------------

% defining expression of Ep and also Tp (if expressed)...
xPR_t     = xPR;
xPR_t(12) = p(1); % ... scaling TX rate of Ep
xPR_t(13) = p(2); % ... scaling TX rate of Tp
xPR_t(10) = p(3); % ... scaling TX rate of E


% --- Calculate objective values -----------------------------------------

% simulate batch culture
pp = m_SimBatchCult_ProdPerf(hPR,xPR_t,SysTopol,y0,1,0,[],[]); % pp = [vProd,pYield]

% output:
obj = (-1)*pp(objID); % this is minimised, which effectively maximises the selected objective objective: 1 = vProd and 2 = pYield
