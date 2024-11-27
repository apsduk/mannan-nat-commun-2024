% ------------------------------------------------------------------------
% RUN_model
% ------------------------------------------------------------------------

%% ===== SET UP ===========================================================
clear('all'); close('all'); addpath(pwd);

% --- Known parameters ----------------------------------------------------
sS0 = 0.5; cultvol = 1.25;

% --- Host export reaction ------------------------------------------------
vX = 726; kX = 1e3;

% --- Time span -----------------------------------------------------------
tmax = 7*24*60;

% --- Unknown bounds for optimisation -------------------------------------
sTX_Z  = [1e-3, 2]; % scaling max TX rate
K_Z    = [1e-6, 1];  % biosensor affinity to bind to inhibit TX of Ep
t_ind  = [1, 24*60];  % range of the time points when inducer is added

% --- Load model parameters -----------------------------------------------
[hPR0, xPR0] = model_params(sS0, vX, kX, cultvol);

% --- Choose solver -------------------------------------------------------
odesolver = @ode15s;

odetolvector = [1e-6, 1e-9, 1e-12];
for i = 1:length(odetolvector)
    odeoptions{i} = odeset('AbsTol', odetolvector(i)*ones(35,1), 'RelTol', odetolvector(i), 'NonNegative', 1:35);
    odeoptions{i}.warning = 0;
end

% --- Make initial conditions ---------------------------------------------
OD600 = 0.05; % inoculated OD
N0  = 1e6;
xS0 = ((10*cultvol)/180)*6.02e23; % amount of glucose in media = 10g/L in 1.25L working vol in 3L vessel
xI0 = 1e23;
Y0scales = [N0 xS0 xI0];

%% ===== BUILD MODEL TO SIMULATE ==========================================

% --- Define topology of interest -----------------------------------------
%          [pre, xp, regT, regE, regEp, regTp, regTF]
SysTopol = [  1,  2,    8,    1,    -1,     0,     0];

% --- Choose parameters ---------------------------------------------------
chooseOpts = [1 0.5 1 0.5 0.1 0.1 300];

% --- Figure set ups ------------------------------------------------------
cgray = [0.8 0.8 0.8];
linewidth = 1;
fontsize = 12;
M0 = 1e8;

% --- Choose solver -------------------------------------------------------
odesolver = @ode15s;

odetolvector = [1e-6];
for i = 1:length(odetolvector)
    odeoptions{i} = odeset('AbsTol', odetolvector(i)*ones(35,1), 'RelTol', odetolvector(i), 'NonNegative', 1:35);
    odeoptions{i}.warning = 0;
end

%% ===== SET UP SIMULATION PROBLEM ========================================
% --- Build bounds --------------------------------------------------------
[lb, ub, xidx, xnames, topodescr] = m_BuildBounds(SysTopol, sTX_Z, K_Z, t_ind);
disp([topodescr])

% --- Print parameters ----------------------------------------------------
for i = 1:length(xnames)
    disp([xnames{i},' = ',num2str(chooseOpts(i))]);
end

% --- Update topo to set con expression -----------------------------------
SysTopol(SysTopol == 8) = 0;
SysTopol(SysTopol == 9) = 0;

% --- Write Jacobian file -------------------------------------------------
jname = writejacobianfile('m_BatchCultModel_DC', ['jacode_temp'], 1, ones(35,1), hPR0, xPR0, Y0scales, SysTopol);

% --- Add name to odeoptions ----------------------------------------------
for i = 1:length(odeoptions)
    odeoptions{i}.jacautogenname = jname;
end

% --- Build model ---------------------------------------------------------
hPR = hPR0;
xPR = xPR0;
xPR(xidx) = chooseOpts(1:end-1);
tind = chooseOpts(end);

% --- Simulate model ------------------------------------------------------
idx = 0;
flag = 0;
while flag == 0
    idx = idx + 1;
    [ProdPerf, lambda, TXrates, TLrates, gammaX, protmass, ribomass, fluxes, summass, T, Y, dY_by_dt, flag] = m_SimBatchCult_SystDynamics(hPR, xPR, Y0scales, SysTopol, tind, tmax, odesolver, odeoptions{idx});
    if idx == length(odeoptions)
        flag = 1;
    end
end

%% ===== PLOT FIGURE ======================================================

% --- Get species out of Y vector -----------------------------------------
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
mEp = Y(:,24); cEp = Y(:,25); pEp = Y(:,26); iP = Y(:,27); % synthesis pathway enzyme
mTp = Y(:,28); cTp = Y(:,29); pTp = Y(:,30);               % product synthesis transporter
mTF = Y(:,31); cTF = Y(:,32); pTF = Y(:,33);               % inducible biosensor for synthetic inducer
iI  = Y(:,34); sC  = Y(:,35);                              % inducer and sequestered complex

% --- Concert time to hours -----------------------------------------------
T = T/60;

% --- Get induction time --------------------------------------------------
tidx = sum(T < tind/60);

%% ===== MAKE FIGURE ======================================================
fplot = figure('Units', 'centimeters', 'Position', [5 5 30 17]);
cmap = lines(8);

% --- Plot substrate and product ------------------------------------------
subplot(1, 2, 2);
yyaxis('left'); hold('on');
plot(T, xS, '-', 'Color', cmap(1,:), 'LineWidth', 2);
plot(T, xP, '-', 'Color', cmap(2,:), 'LineWidth', 2);
set(gca, 'YColor', 'k');
ylabel('External env. (molecules)');
yyaxis('right'); hold('on');
plot(T, N, '-', 'Color', 'k', 'LineWidth', 2);
set(gca, 'YColor', 'k');
ylabel('N (cells)');
xlim([min(T), max(T)]);
xlabel('Time (min)');
plot(T(tidx), min(ylim), '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
plot(T(tidx), max(ylim), 'vk', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
plot([T(tidx); T(tidx)], [min(ylim); max(ylim)], ':k', 'LineWidth', linewidth);
legend('Ext. S','Ext. P', 'N', 'Orientation', 'horizontal', 'Location', 'northoutside');
% title(topodescr);
set(gca, 'Box', 'on', 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', 1);

% --- Metabolic rates -----------------------------------------------------
subplot(1, 2, 1); 
yyaxis('left'); hold('on');
plot(T, fluxes(:,2)./fluxes(1,2), '-', 'Color', cmap(4,:),  'LineWidth', 2);
plot(T, sum(fluxes(:,[3,4]),2)./sum(fluxes(1,[3,4]),2), '-', 'Color', cmap(2,:),  'LineWidth', 2);
set(gca, 'YColor', 'k');
xlim([min(T), max(T)]);
set(gca, 'Box', 'on', 'FontSize', fontsize, 'FontWeight', 'bold', 'LineWidth', 1);
ylabel('Metabolic rates (molecules per min)');
yyaxis('right');
plot(T, lambda, '-k', 'LineWidth', 2);
ylabel('Growth rate (per min)');
set(gca, 'YColor', 'k');
plot(T(tidx), min(ylim), '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
plot(T(tidx), max(ylim), 'vk', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
plot([T(tidx); T(tidx)], [min(ylim); max(ylim)], ':k', 'LineWidth', linewidth);
if SysTopol(1) == 1
    legend('S --> "energy"', 'S --> P', 'Growth rate', 'Orientation', 'horizontal', 'Location', 'northoutside');
elseif SysTopol(1) == 2
    legend('S --> "energy"', 'S --> P', 'Growth rate', 'Orientation', 'horizontal', 'Location', 'northoutside');
end
subplot(1, 2, 1); set(gca, 'Box', 'on', 'PlotBoxAspectRatio', [1 1 1]); xticks([0:6:max(T)]); yyaxis('left'); yyaxis('right'); 
subplot(1, 2, 2); set(gca, 'Box', 'on', 'PlotBoxAspectRatio', [1 1 1]); xticks([0:6:max(T)]); yyaxis('left'); yyaxis('right');
xlabel('Time (min)');
