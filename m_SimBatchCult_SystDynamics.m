% -------------------------------------------------------------------------
% SIMULATE BATCH CULTURE OF SYSTEM WITH DYNAMIC CONTROLLER AND MEASURE
% PRODUCTION PERFORMANCE
% -------------------------------------------------------------------------
% Alexander P.S. Darlington
% -------------------------------------------------------------------------
function [ProdPerf, lambda, TXrates, TLrates, gammaX, protmass, ribomass, fluxes, summass, T, Y, dY_by_dt, flag] = m_SimBatchCult_SystDynamics(hPR, xPR, y0scales, SysTopol, tind, tmax, odesolver, odeoptions)

%% ===== SIMULATE MODEL ===================================================
[ProdPerf, T, Y, flag, tsol, ysol] = m_SimBatchCult_ProdPerf(hPR, xPR, y0scales, SysTopol, tind, tmax, odesolver, odeoptions);

%% ===== ITERATE OVER T AND Y =============================================
% --- Make space ----------------------------------------------------------
dY_by_dt = zeros(length(T), 35);
lambda   = zeros(length(T),  1);
TXrates  = zeros(length(T),  9);
TLrates  = zeros(length(T),  8);
gammaX   = zeros(length(T),  1);
protmass = zeros(length(T),  8);
summass  = zeros(length(T),  1);
ribomass = zeros(length(T),  1);
fluxes   = zeros(length(T),  6);

% --- Iterate over T and Y ------------------------------------------------
for t = 1:length(T)
    
    % --- Simulate dY at point t ------------------------------------------
    [dY, lambda(t), TXrates(t,:), TLrates(t,:), gammaX(t,:), protmass(t,:), summass(t,:), ribomass(t,:), fluxes(t,:)] = m_BatchCultModel_DC(T(t), Y(t,:), hPR, xPR, y0scales, SysTopol);
    dY_by_dt(t,:) = dY';
    
end

% --- Get outputs of interest ---------------------------------------------
N   = Y(:, 1); % total biomass
xS  = Y(:, 2); % substrate in media
xP  = Y(:, 3); % product in media
xI  = Y(:, 4); % inducer conc in media
iS  = Y(:, 5); ee  = Y(:, 6);                % internal substrate (is) and energY (ee)
mT  = Y(:, 7); cT  = Y(:, 8); pT  = Y(:, 9); % transporter
mE  = Y(:,10); cE  = Y(:,11); pE  = Y(:,12); % enzymes
mH  = Y(:,13); cH  = Y(:,14); pH  = Y(:,15); % approx fixed amount of host proteins (q-fraction)
mX  = Y(:,16); cX  = Y(:,17); pX  = Y(:,18); % export reaction
mR  = Y(:,19); cR  = Y(:,20); xR  = Y(:,21); rr = Y(:,22); pR = Y(:,23); % ribosomes
mEp = Y(:,24); cEp = Y(:,25); pEp = Y(:,26); iP = Y(:,27);               % synthesis pathway enzyme
mTp = Y(:,28); cTp = Y(:,29); pTp = Y(:,30);                             % product synthesis transporter
mTF = Y(:,31); cTF = Y(:,32); pTF = Y(:,33);                             % inducible biosensor for synthetic inducer
iI  = Y(:,34); sC  = Y(:,35);                                            % inducer and sequestered complex


end

