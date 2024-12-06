% ------------------------------------------------------------------------
% OBJECTIVES OF MULTIOBJECTIVE OPTIMIZATION PROBLEM 1
%
% Objectives to be minimized are:
%   - obj1 = -1 * growth rate
%   - obj2 = -1 * product synthesis flux
% ------------------------------------------------------------------------

function obj = m_MultiObjs_1(p,hPR,xPR,y0,SysTopol)


% --- Redefine parameters ------------------------------------------------

% defining expression of Ep and also Tp (if expressed)...
xPR_t     = xPR;
xPR_t(12) = p(1); % ... scaling TX rate of Ep
xPR_t(13) = p(2); % ... scaling TX rate of Tp
xPR_t(10) = p(3); % ... scaling TX rate of E
xPR_t(14) = p(4); % ... scaling TX rate of TF


% --- Calculate objective values -----------------------------------------

% calculate steady state of the system for given param set:
[~,growthr,fluxes] = m_steadystate(hPR,xPR_t,y0,SysTopol);

% objective 1:
obj1 = (-1) * growthr;

% objective 2:
% obj2 = (-1) * (fluxes(3) + fluxes(4));
obj2 = (-1) * fluxes(5);

% output:
obj = [obj1,obj2];