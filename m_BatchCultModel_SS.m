% ------------------------------------------------------------------------
% GROWTH AND PRODUCTION IN BATCH CULTURE
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% Edits Alexander P.S. Darlington
% ------------------------------------------------------------------------

function [dY, lambda, TXrates, TLrates, gammaX, protmass, summass, ribomassfrac, fluxes] = m_BatchCultModel_SS(T, Y, hPR, xPR, y0scales, SysTopol)

% --- Calculate derivative ------------------------------------------------
[dY, lambda, TXrates, TLrates, gammaX, protmass, summass, ribomassfrac, fluxes] = m_BatchCultModel_DC(T, Y, hPR, xPR, y0scales, SysTopol);

% --- Make extracellular reactions zero -----------------------------------
dY(1:4) = zeros(4,1);


end