% ------------------------------------------------------------------------
% ANALYSIS 2 - MULTIOBJECTIVE OPTIMIZATION TO MAX 
%              (i)  GROWTH AND PRODUCT SYNTHESIS FLUX
%              (ii) VOLUMETRIC PRODUCTIVITY AND YIELD
% ------------------------------------------------------------------------

function A1_MultiobjOpt(is,xporter,TuneEpAndE,fig)


figure(fig); clf
for i = 1:6
    figure(fig+i); clf
end

% ------------------------------------------------------------------------
% [A] SOLVE PROBLEM (i)
% ------------------------------------------------------------------------

% (a) Define params, bounds on those to explore over, and initial conds
% ------------------------------------------------------------------------

% Define nominal parameters:
[n_hPR,n_xPR] = model_params;

% Define prod precursor and topology TF-mediated dynamic control circuit:
%           is?  X?       regT regE regX regEp regTp regTF
SysTopol = [is,  xporter, 0,   0,   0,   0,    0,    0    ]; % see m_defSysAndTopols.m for details of vector definition

% Ranges of ...
% ... scaling to max TX rate of Ep:
sTX_Ep = [1e-3,1];
% ... scaling to max TX rate of Tp:
if xporter == 1
    sTX_Tp = [0,0];
elseif xporter == 2
    sTX_Tp = [1e-3,1];
end
% ... scale to TX rate of host enzyme E:
if TuneEpAndE == 1
    sTX_E = [1e-3,1];
else
    sTX_E = [1,1];
end
% ... scale to TX rate of TF: ############################################
% sTX_TF = [1e-3,1];
sTX_TF = [0,0];


% Defining lower and upper bounds to the two params to explore over:
lb = [sTX_Ep(1),sTX_Tp(1),sTX_E(1),sTX_TF(1)];
ub = [sTX_Ep(2),sTX_Tp(2),sTX_E(2),sTX_TF(2)];

% Define initial conditions, i.e. steady state w/o synthesis p/w:
x0 = m_ssFor0pwExp(n_hPR,n_xPR,SysTopol);


% (b) Run the multiobjective optimization problem
% ------------------------------------------------------------------------

% Solve optimization problem:
options = optimoptions(@gamultiobj,'Display','iter','PopulationSize',200,'MaxStallGenerations',10,'UseParallel',true,'MaxGenerations',100,'PlotFcn',@gaplotpareto);
[pOpts_t,objVals_t,exitflag] = gamultiobj(@(p)m_MultiObjs_1(p,n_hPR,n_xPR,x0,SysTopol),length(lb),[],[],[],[],lb,ub,options);
objVals_t = (-1) * objVals_t;
disp(exitflag)

% Reorder objVals in ascending order of growth rate:
[~,idx] = sort(objVals_t(:,1),'ascend');

% Redefine objective vals and opt params in ascending order of growth:
objVals = objVals_t(idx,:);
pOpts = pOpts_t(idx,:);


% (c) Calculate vol. productivity and yield for each optimal design
% ------------------------------------------------------------------------

ProdVsYld = zeros(size(objVals));
for i = 1:length(objVals(:,1))
    [~,ov] = m_MultiObjs_2(pOpts(i,:),n_hPR,n_xPR,x0,SysTopol,[1,1]);
    
    ProdVsYld(i,:) = ov;
end

save('MultiObj_1', 'pOpts', 'objVals', 'ProdVsYld');

%%
% (d) Plotting the Pareto solutions
% ------------------------------------------------------------------------

% Plot of Pareto fronts:
% ... of synthesis vs growth rates:
figure(fig); clf
subplot(2,1,1)
% plot(objVals(:,1),objVals(:,2),'ko-','MarkerSize',8)
scatter(objVals(:,1),objVals(:,2),45,objVals(:,2),'x','LineWidth',2)
xlabel('Obj_1 | Growth rate (1/min)')
ylabel('Obj_2 | Synthesis flux r_{T_p} (mol/cell/min)')
colormap('summer')
c = colorbar;
c.Label.String = 'Synthesis flux';
c.Label.FontSize = 10;
title('Results from maximising synthesis and growth')
% ... of productivity vs yield:
axis([0 0.04 0 2.5e7])

subplot(2,1,2)
% plot(ProdVsYld(:,2),ProdVsYld(:,1),'ko','MarkerSize',8)
scatter(ProdVsYld(:,2),ProdVsYld(:,1),45,objVals(:,2),'x','LineWidth',2)
xlabel('Obj_1 | Yield (g-prod/g-subs)')
ylabel('Obj_2 | Vol. Productivity (molecs/min')
c = colorbar;
c.Label.String = 'Synthesis flux';
c.Label.FontSize = 10;
axis([0 1 0 4e19])

% Select and mark the three example designs Plot of time-course sims of three designs of different performances:
idx_des = round([0.20,0.55,0.91]*length(objVals(:,1)));

% figure(fig)
subplot(2,1,1)
hold on
plot(objVals(idx_des(1),1),objVals(idx_des(1),2),'ko','MarkerSize',14)
plot(objVals(idx_des(2),1),objVals(idx_des(2),2),'o','MarkerSize',14,'Color',[0.5,0.5,0.5])
plot(objVals(idx_des(3),1),objVals(idx_des(3),2),'o','MarkerSize',14,'Color',[0.7,0.7,0.7])
hold off
subplot(2,1,2)
hold on
plot(ProdVsYld(idx_des(1),2),ProdVsYld(idx_des(1),1),'ko','MarkerSize',14)
plot(ProdVsYld(idx_des(2),2),ProdVsYld(idx_des(2),1),'o','MarkerSize',14,'Color',[0.5,0.5,0.5])
plot(ProdVsYld(idx_des(3),2),ProdVsYld(idx_des(3),1),'o','MarkerSize',14,'Color',[0.7,0.7,0.7])
hold off


% Plot of parameter values:

figure(fig+1); clf
subplot(3,1,1)
% plot(pOpts(:,1),objVals(:,2),'ko','MarkerSize',6,'MarkerFaceColor',[0.5,0.5,0.5])
scatter(pOpts(:,1),objVals(:,2),45,objVals(:,2),'x','LineWidth',2)
ylabel('Synthesis flux (mol/cell/min)')
xlabel('Scale TX rate of E_p')
xlim(sTX_Ep); set(gca,'XScale','log')
colormap('summer')

subplot(3,1,2)
% plot(pOpts(:,2),objVals(:,2),'ko','MarkerSize',6,'MarkerFaceColor',[0.5,0.5,0.5])
scatter(pOpts(:,2),objVals(:,2),45,objVals(:,2),'x','LineWidth',2)
ylabel('Synthesis flux (mol/cell/min)')
xlabel('Scale TX rate of T_p')
xlim(sTX_Tp); set(gca,'XScale','log')

subplot(3,1,3)
% plot(pOpts(:,3),objVals(:,2),'ko','MarkerSize',6,'MarkerFaceColor',[0.5,0.5,0.5])
scatter(pOpts(:,3),objVals(:,2),45,objVals(:,2),'x','LineWidth',2)
ylabel('Synthesis flux (mol/cell/min)')
xlabel('Scale TX rate of E')
set(gca,'YScale','log')
xlim(sTX_E); set(gca,'XScale','log')


figure(fig+2); clf
subplot(3,1,1)
scatter(pOpts(:,1),pOpts(:,2),45,objVals(:,2),'x','LineWidth',2)
xlabel('Scale TX rate of E_p')
ylabel('Scale TX rate of T_p')
xlim(sTX_Ep); set(gca,'XScale','log')
ylim(sTX_Tp); set(gca,'YScale','log')
colormap('summer')

subplot(3,1,2)
scatter(pOpts(:,2),pOpts(:,3),45,objVals(:,2),'x','LineWidth',2)
xlabel('Scale TX rate of T_p')
ylabel('Scale TX rate of E')
xlim(sTX_Tp); set(gca,'XScale','log')
ylim(sTX_E); set(gca,'YScale','log')

subplot(3,1,3)
scatter(pOpts(:,3),pOpts(:,1),45,objVals(:,2),'x','LineWidth',2)
xlabel('Scale TX rate of E')
ylabel('Scale TX rate of E_p')
xlim(sTX_E); set(gca,'XScale','log')
ylim(sTX_Ep); set(gca,'YScale','log')


% Plot of time-course sims of three designs of different performances:
figure(fig+3); clf
for i = 1:length(idx_des)
    % simulate culture time-course for a defined strain design:
    n_xPR_t     = n_xPR;  % defining design
    row_pOpt    = idx_des(i);
    n_xPR_t(12) = pOpts(row_pOpt,1); % ... scaling TX rate of Ep
    n_xPR_t(13) = pOpts(row_pOpt,2); % ... scaling TX rate of Tp
    n_xPR_t(10) = pOpts(row_pOpt,3); % ... scaling TX rate of E
    n_xPR_t(14) = pOpts(row_pOpt,4); % ... scaling TX rate of TF
    [~,simT,simY] = m_SimBatchCult_ProdPerf(n_hPR,n_xPR_t,SysTopol,x0,1,0,[],[]);

    % plotting time-course:
    subplot(length(idx_des),1,i)
    ax = gca();

    yyaxis left
    plot(simT,simY(:,2),'-','LineWidth',1,'Color',[0.5,0.5,0.5])
    hold on
    plot(simT,simY(:,3),'-','LineWidth',5,'Color',[0.4660 0.6740 0.1880])
    hold off
    ylabel('Number of molecules')
    set(gca,'YScale','lin'); ylim([0,4.5e22])
    
    yyaxis right
    plot(simT,simY(:,1),'k-','LineWidth',2)
    ylabel('Number of cells')
    set(gca,'YScale','log'); ylim([1e6,1e16])

    ax.YAxis(1).Color = [0.4660 0.6740 0.1880];
    ax.YAxis(2).Color = [0,0,0];
    xlim([0,2500])
    % legend({'media substrate','product','biomass'},'Location','north')
    xlabel('Time (mins)')
end


%% ------------------------------------------------------------------------
% [B] SOLVE PROBLEM (ii)
% ------------------------------------------------------------------------

% (a) Define params, bounds on those to explore over, and initial conds
% ------------------------------------------------------------------------

% Define nominal parameters:
[n_hPR,n_xPR] = model_params;

% Define prod precursor and topology TF-mediated dynamic control circuit:
%           is?  X?       regT regE regX regEp regTp regTF
SysTopol = [is,  xporter, 0,   0,   0,   0,    0,    0    ]; % see m_defSysAndTopols.m for details of vector definition

% Ranges of ...
% ... scaling to max TX rate of Ep:
sTX_Ep = [1e-3,1e0];
% ... scaling to max TX rate of Tp:
if xporter == 1
    sTX_Tp = [0,0];
elseif xporter == 2
    sTX_Tp = [1e-3,1e0];
end
% ... scale to TX rate of host enzyme E:
if TuneEpAndE == 1
    sTX_E = [1e-3,1e0];
else
    sTX_E = [1,1];
end
% ... scale to TX rate of TF: ############################################
% sTX_TF = [1e-3,1];
sTX_TF = [0,0];

% Defining lower and upper bounds to the two params to explore over:
lb = [sTX_Ep(1),sTX_Tp(1),sTX_E(1),sTX_TF(1)];
ub = [sTX_Ep(2),sTX_Tp(2),sTX_E(2),sTX_TF(2)];

% Define initial conditions, i.e. steady state w/o synthesis p/w:
x0 = m_ssFor0pwExp(n_hPR,n_xPR,SysTopol);


% (b) Run the multiobjective optimization problem - max prod and yield
% ------------------------------------------------------------------------

% Solve optimization problem:
options = optimoptions(@gamultiobj,'Display','iter','PopulationSize',200,'MaxStallGenerations',10,'UseParallel',true,'MaxGenerations',100,'PlotFcn',@gaplotpareto);
[pOpts_t,ProdYields_t,exitflag] = gamultiobj(@(p)m_MultiObjs_2(p,n_hPR,n_xPR,x0,SysTopol,[1,1]),length(lb),[],[],[],[],lb,ub,options);
ProdYields_t = (-1) * ProdYields_t;
disp(exitflag)

% Reorder objVals in ascending order of productivity:
[~,idx] = sort(ProdYields_t(:,1),'ascend');

% Redefine objective vals and opt params in ascending order of growth:
ProdYields = ProdYields_t(idx,:);
pOpts = pOpts_t(idx,:);


% (c) Calculate synthesis and growth rates for each optimal design
% ------------------------------------------------------------------------

SynthGrowth = zeros(size(ProdYields));
for i = 1:length(ProdYields(:,1))
    objs_t = m_MultiObjs_1(pOpts(i,:),n_hPR,n_xPR,x0,SysTopol);
    
    SynthGrowth(i,:) = (-1) * objs_t;
end

save('MultiObj_2', 'pOpts', 'ProdYields', 'SynthGrowth');



% (d) Plotting the Pareto solutions
% ------------------------------------------------------------------------

% Plot of Pareto fronts:
% ... of synthesis vs growth rates:
figure(fig+4); clf
subplot(2,1,2)
scatter(ProdYields(:,2),ProdYields(:,1),50,ProdYields(:,2),'filled')
xlabel('Obj_1 | Yield (g-prod/g-subs)')
ylabel('Obj_2 | Vol. Productivity (molecs/min)')
colormap('autumn')
c = colorbar;
c.Label.String = 'Yield (g-prod/g-subs)';
c.Label.FontSize = 10;
axis([0 1 0 4e19])

% ... of productivity vs yield:
subplot(2,1,1)
scatter(SynthGrowth(:,1),SynthGrowth(:,2),50,ProdYields(:,2),'filled')
xlabel('Obj_1 | Growth rate (1/min)')
ylabel('Obj_2 | Synthesis flux r_{T_p} (mol/cell/min)')
c = colorbar;
c.Label.String = 'Yield (g-prod/g-subs)';
c.Label.FontSize = 10;
title('Results from maximising vol. prod. and yield')
axis([0 0.04 0 2.5e7])


% Plot of parameter values:

figure(fig+5); clf
subplot(3,1,1)
scatter(pOpts(:,1),ProdYields(:,2),50,ProdYields(:,2),'filled')
ylabel('Yield (g-prod/g-subs)')
xlabel('Scale TX rate of E_p')
xlim(sTX_Ep); set(gca,'XScale','log')
ylim([0,1])
colormap('autumn')

subplot(3,1,2)
scatter(pOpts(:,2),ProdYields(:,2),50,ProdYields(:,2),'filled')
ylabel('Yield (g-prod/g-subs)')
xlabel('Scale TX rate of T_p')
xlim(sTX_Tp); set(gca,'XScale','log')
ylim([0,1])

subplot(3,1,3)
scatter(pOpts(:,3),ProdYields(:,2),50,ProdYields(:,2),'filled')
ylabel('Yield (g-prod/g-subs)')
xlabel('Scale TX rate of E')
set(gca,'YScale','log')
xlim(sTX_E); set(gca,'XScale','log')
ylim([0,1])


figure(fig+6); clf
subplot(3,1,1)
scatter(pOpts(:,1),pOpts(:,2),50,ProdYields(:,2),'filled')
xlabel('Scale TX rate of E_p')
ylabel('Scale TX rate of T_p')
xlim(sTX_Ep); set(gca,'XScale','log')
ylim(sTX_Tp); set(gca,'YScale','log')
colormap('autumn')

subplot(3,1,2)
scatter(pOpts(:,2),pOpts(:,3),50,ProdYields(:,2),'filled')
xlabel('Scale TX rate of T_p')
ylabel('Scale TX rate of E')
xlim(sTX_Tp); set(gca,'XScale','log')
ylim(sTX_E); set(gca,'YScale','log')

subplot(3,1,3)
scatter(pOpts(:,3),pOpts(:,1),50,ProdYields(:,2),'filled')
xlabel('Scale TX rate of E')
ylabel('Scale TX rate of E_p')
set(gca,'YScale','log')
xlim(sTX_E); set(gca,'XScale','log')
ylim(sTX_Ep); set(gca,'YScale','log')
