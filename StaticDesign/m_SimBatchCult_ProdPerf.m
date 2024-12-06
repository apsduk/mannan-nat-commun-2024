% ------------------------------------------------------------------------
% SIMULATE BATCH CULTURE OF SYSTEM AND MEASURING PRODUCTION PERFORMANCE
% ------------------------------------------------------------------------

function [ProdPerf,T,Y,ss0,GR] = m_SimBatchCult_ProdPerf(hPR,xPR,SysTopol,ss0,t_ind,Ix_ind,MeasGR,CultSimFig)


% --- Set up for simulation ----------------------------------------------

OD600 = 0.05; % inoculated OD
mLvol = 1250; % working volume of 1.25L
% X0  = ((OD600 - 0.01)/0.446)*1e9*(mLvol); %9e7; % inoculating 1.25L working vol to 0.05 OD600 of biomass, based on Yap2019
X0  = 1e6; % ~ OD600 = 0.010000357
xS0 = ((10*1.25)/180)*6.02e23; % amount of glucose in media = 10g/L in 1.25L working vol in 3L vessel
xP0 = 0; % product in media
xI0 = 0; % initial amount of inducer added
SS0 = [X0,xS0,xP0,xI0,ss0];


% --- Simulate batch culture ---------------------------------------------

% Define event func. for ODE to terminate sims once substrate depletes:
opts    = odeset('Events',@m_Event_SubsRunOut);

% Stage 1 - culture w/o inducer:
tspan1  = [0 t_ind]; % almost 7 days time length
ic1     = SS0;
[t1,y1] = ode15s(@(t,y)m_BatchCultModel(t,y,hPR,xPR,SysTopol),tspan1,ic1,opts);

% Stage 2 - addition of inducer:
if y1(end,2) >= 0 % if there is substrate left over ...
    tspan2  = t1(end) + [0 60*24*7];
    ic2     = y1(end,:);
    ic2(4)  = Ix_ind; % 5e23; % amount of inducer added at time t_ind (in molecules)
    [t2,y2] = ode15s(@(t,y)m_BatchCultModel(t,y,hPR,xPR,SysTopol),tspan2,ic2,opts);
    
    switchInd = 1; % indicates if switch induction happened
else
    t2 = [];
    y2 = [];
    switchInd = 0; % indicates if switch induction happened
end

% Full time courses:
T = [t1; t2];
Y = [y1; y2];

% Measure growth rate at each time point:
if ~isempty(MeasGR)
    GR = zeros(length(T),1);
    for i = 1:length(T)
        [~,gr] = m_BatchCultModel(T(i),Y(i,:),hPR,xPR,SysTopol);
        GR(i) = gr;
    end
end


% --- Plotting -----------------------------------------------------------

if ~isempty(CultSimFig)
    figure(CultSimFig); clf
    
    subplot(3,2,1) % time course of substrate and inducer
    [aX,L1,L2] = plotyy(T,Y(:,2),T,Y(:,4));
    xlabel('Time (mins)');
    ylabel(aX(1),{'Substrate';'(molecules)'})
    ylabel(aX(2),{'Inducer (ext)';'(molecules)'})
    L1.LineWidth = 2; L2.LineWidth = 2;
    xlim(aX, [min(T) max(T)])
    
    subplot(3,2,3) % time course of population and product
    [aX,L1,L2] = plotyy(T,Y(:,1),T,Y(:,3));
    xlabel('Time (mins)')
    ylabel(aX(1),{'Population';'(num. of cells)'})
    ylabel(aX(2),{'Product';'(molecules)'})
    L1.LineWidth = 2; L2.LineWidth = 2;
    set(aX(1),'YScale','log','YTickMode','auto')
    xlim(aX, [min(T) max(T)])
    
    subplot(3,2,5) % time course of product
    plot(T,Y(:,3),'-','LineWidth',2)
    xlabel('Time (mins)'); ylabel({'Product';'(molecules)'})
%     set(gca,'YScale','log')
    ylim([0.001,1]*max(Y(:,3)))
    xlim([min(T) max(T)])
    
    subplot(3,2,2) % time course of inducer in cell
    plot(T,Y(:,31),'-','LineWidth',2)
    xlabel('Time (mins)'); ylabel({'Inducer (int)';'(molecules/cell)'})
    xlim([min(T) max(T)])
    
    subplot(3,2,4) % time course of sequestered complex
    plot(T,Y(:,32),'-','LineWidth',2)
    xlabel('Time (mins)'); ylabel({'Seq.Complex';'(molecules/cell)'})
    xlim([min(T) max(T)])
    
    subplot(3,2,6) % time course of growth rate
    plot(T,GR,'-','LineWidth',2)
    xlabel('Time (mins)'); ylabel({'Growth rate';'(1/min)'})
    xlim([min(T) max(T)])
end


% --- Measuring production performance -----------------------------------

% if switchInd == 1
    % Volumetric productivity (total prod / total batch culture time):
    vProd  = Y(end,3)/T(end);
    
    % Product yield (total prod / total substrate consumed):
    pYield = Y(end,3)/(Y(1,2)-Y(end,2));
    
    ProdPerf = [vProd,pYield];
% else
%     ProdPerf = [0,0];
% end
