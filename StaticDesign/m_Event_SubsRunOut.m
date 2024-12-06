% ------------------------------------------------------------------------
% EVENT LOCATOR THAT STOPS BATCH CULTURE SIM WHEN SUBSTRATE RUNS OUT
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% ------------------------------------------------------------------------


function [val,isterminal,dir] = m_Event_SubsRunOut(t,y) %#ok<INUSL>

val = y(1) - (1/1); % i.e. looking at concentration of 1nM substrate in 1.25L working volume
isterminal = 1;
dir = 0;