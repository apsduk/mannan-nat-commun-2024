% ------------------------------------------------------------------------
% CARRY OUT PRODUCT YIELD OPTIMISATION
% ------------------------------------------------------------------------

%% ===== SET UP ===========================================================
clear('all'); close('all');

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
[hPR, xPR] = model_params(sS0, vX, kX, cultvol);

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

% --- Options for optimisation --------------------------------------------
initnpop =  100;
ganpop   =  500;
nga      =  500;
nmultiga = 1000;

%% ===== SET UP OPTIMISATIONS =============================================

% --- Single Obj GA options -----------------------------------------------
gaoptions = optimoptions(@ga, 'UseParallel',true, ...
    'PopulationSize', initnpop, 'MaxGenerations', nga, 'MaxStallGenerations', ceil(nga/10));

% --- Multi GA options ----------------------------------------------------
multigaoptions = optimoptions(@gamultiobj, 'UseParallel', true, ...
    'PopulationSize', ganpop, 'MaxGenerations', nmultiga, 'MaxStallGenerations', ceil(nmultiga/10));

%% ===== SET TOPOLOGIES ===================================================
allSysTopols = m_defSysAndTopols;

disp(['There are ',num2str(length(allSysTopols(:,1))),' topologies to test']);

%% ===== MAKE RESULTS FOLDER ==============================================

% --- Make folders --------------------------------------------------------
addpath(pwd);
fname = 'example_optimisation_results';
mkdir(fname);
cd(fname);
fname = pwd;

pwd

% --- Initiate random number seed -----------------------------------------
rngsettings = rng('shuffle');

% --- Save results --------------------------------------------------------
save('simulation_setup.mat')

% --- Open parallel pool --------------------------------------------------
if isunix
    parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')))
end

%% ===== SELECT TOPOLOGY INDEX OF INTEREST ================================
topoid = 56;

%% ===== OPTIMISE FRONT ===================================================
pwd

% --- Get topology and set up folders ---------------------------------
disp(['TopoID ',num2str(topoid)]);
addpath(pwd);
SysTopol = allSysTopols(topoid,:);

% --- Build bounds ----------------------------------------------------
[lb, ub, xidx, xnames, topodescr] = m_BuildBounds(SysTopol, sTX_Z, K_Z, t_ind);
disp(['TopoID ',num2str(topoid),' | ',topodescr])

% --- Update topo to set con expression -------------------------------
SysTopol(SysTopol == 8) = 0;
SysTopol(SysTopol == 9) = 0;

% --- Write Jacobian file ---------------------------------------------
jname = writejacobianfile('m_BatchCultModel_DC', ['jacode_topoid_',num2str(topoid)], 1, ones(35,1), hPR, xPR, Y0scales, SysTopol);

% --- Add name to odeoptions ------------------------------------------
for i = 1:length(odeoptions)
    odeoptions{i}.jacautogenname = jname;
end

% --- Build volumetric productivity objective function ----------------
vprodObjfun = @(p) m_VolProdObj_v2(p, hPR, xPR, Y0scales, SysTopol, tmax, xidx, odesolver, odeoptions);

% --- Build volumetric productivity objective function ----------------
pyieldObjfun = @(p) m_ProdYieldObj_v2(p, hPR, xPR, Y0scales, SysTopol, tmax, xidx, odesolver, odeoptions);

% --- Carryout optimisation for prod yield ----------------------------
disp(['Carry out optimisation for pYield using ga | ',datestr(now(),'yyyy-mm-dd-HHMM')]);
[pOpts, fval, exitflag, outputstruct, population, scores] = ga(pyieldObjfun, length(lb), [], [], [], [], lb, ub, [], [], gaoptions);
disp(exitflag)
gaPYieldStruct.x = pOpts;
gaPYieldStruct.fval = fval;
gaPYieldStruct.exitflag = exitflag;
gaPYieldStruct.output = outputstruct;
gaPYieldStruct.population = population;
gaPYieldStruct.scores = scores;
maxpyield = abs(fval);

% --- Carryout optimisation for vol prod ------------------------------
disp(['Carry out optimisation for VProd using ga | ',datestr(now(),'yyyy-mm-dd-HHMM')]);
[pOpts, fval, exitflag, outputstruct, population, scores] = ga(vprodObjfun, length(lb), [], [], [], [], lb, ub, [], [], gaoptions);
disp(exitflag)
gaVProdStruct.x = pOpts;
gaVProdStruct.fval = fval;
gaVProdStruct.exitflag = exitflag;
gaVProdStruct.output = outputstruct;
gaVProdStruct.population = population;
gaVProdStruct.scores = scores;
maxvprod = abs(fval);

% --- Calculate scaling factor ----------------------------------------
sfac = [maxvprod, maxpyield];

% --- Build initial population range ----------------------------------
initpOpts = sort([gaVProdStruct.x; gaPYieldStruct.x],1);

% --- Update gaoptions ------------------------------------------------
multigaoptions.InitialPopulationRange = initpOpts;

% --- Build trade-off objective function ------------------------------
frontObjfun = @(p) m_MultiObj_v2(p, hPR, xPR, Y0scales, SysTopol, tmax, sfac, xidx, odesolver, odeoptions);

% --- Carryout multiGA optimisation -----------------------------------
disp(['Carry out optimisation for VProd and Yield using gamultiobj | ',datestr(now(),'yyyy-mm-dd-HHMM')]);
[pOpts, fval, exitflag, outputstruct, population, scores] = gamultiobj(frontObjfun, length(lb), [], [], [], [], lb, ub, multigaoptions);
disp(exitflag)
multigaFrontStruct.x = pOpts;
multigaFrontStruct.fval = fval;
multigaFrontStruct.exitflag = exitflag;
multigaFrontStruct.output = outputstruct;
multigaFrontStruct.population = population;
multigaFrontStruct.scores = scores;

% --- Calculate front points ------------------------------------------
objVals = zeros(length(pOpts(:,1)),2);
for i = 1:length(objVals(:,1))

    % --- Simulate point on front -------------------------------------
    [~, objVals(i,:), T, Y] = frontObjfun(pOpts(i,:));

end

% --- Save and return -------------------------------------------------
save(['topo_',num2str(topoid),'_results.mat']);

% --- Delete parpool ------------------------------------------------------
delete(gcp);

% --- Return --------------------------------------------------------------
cd('..');
