% ------------------------------------------------------------------------
% Steady state of system in absence of expressing heterologous pathway
% ------------------------------------------------------------------------

function ss_0exp = m_ssFor0pwExp(h_params,x_params,SysTopols)

% Define initial conditions:
initConds = zeros(1,31);
initConds([1,2,14,15,16]) = [1,1,1,1,1];

% Calculate approximate steady state concentrations:
ss_0exp = m_steadystate(h_params,x_params,initConds,SysTopols);
