% -------------------------------------------------------------------------
% SIMULATE BATCH CULTURE OF SYSTEM WITH DYNAMIC CONTROLLER AND MEASURE
% PRODUCTION PERFORMANCE
% -------------------------------------------------------------------------
% Ahmad A Mannan
% Edits Alexander P.S. Darlington
% -------------------------------------------------------------------------
function [ProdPerf, T, Y, flag, tsol, ysol] = m_SimBatchCult_ProdPerf(hPR, xPR, y0scales, SysTopol, tind, tmax, odesolver, odeoptions)

% --- Turn warnings on/off ------------------------------------------------
if odeoptions.warning == 1; warning('on'); else; warning('off'); end

% --- Build jacobian ------------------------------------------------------
jacode = str2func(odeoptions.jacautogenname);

% --- Jacobian function ---------------------------------------------------
jacfun = @(T,Y) jacode(T, Y, hPR, xPR, y0scales);

% --- Update ode options --------------------------------------------------
odeoptions.Jacobian = jacfun;

%% ===== SIMULATE INITIAL CONDITIONS ======================================
% --- Run in time for initial conditions generation -----------------------
runintmax = 1e6;

% --- Steady state options ------------------------------------------------
% ssodeoptions = odeset('NonNegative', 1:35);

% --- Create and simulate initial conditions ------------------------------
Y0 = zeros(35,1);
Y0( 1) = 0; % N0
Y0( 2) = 1e4; % /y0scales(2); % xS0
Y0( 5) = 1e6; % /y0scales(2); % iS0
Y0( 6) = 1e6; % /y0scales(2); % ee0
Y0( 9) = 1e2; % pT0
Y0(12) = 1e2; % pE0
Y0(23) = 1e2; % pR0
Y0(33) = 1; % pTF0
[T0, Y0] = odesolver( @(T, Y) m_BatchCultModel_SS(T, Y, hPR, xPR, y0scales, SysTopol), [0, runintmax], Y0, odeoptions);

% --- Update initial conditions -------------------------------------------
Y0 = Y0(end,:)'; Y0(Y0 < 0) = 0;
Y0(1) = y0scales(1); % 1; % N0
Y0(2) = y0scales(2); % 1; % xS0
Y0(3) = 0; % xP0
Y0(4) = 0; % xI0

%% ===== SIMULATE BATCH CUTLURE :: Pre tind ===============================

% --- Simulate pre induction ----------------------------------------------
[t1, y1] = odesolver( @(T, Y) m_BatchCultModel_DC(T, Y, hPR, xPR, y0scales, SysTopol), [0, tind], Y0, odeoptions);

%% ===== SIMULATE BATCH CULTURE :: Post tind ==============================

% --- Update initial conditions -------------------------------------------
Y0 = y1(end,:); Y0(Y0 < 0) = 0;
Y0(4) = y0scales(3); % 1; % add inducer xI0

% --- Check if we can induce the system -----------------------------------
if y1(end,2) > 0 % then substrate is left so can induce
    
    % --- Simulate post induction -----------------------------------------
    [t2, y2] = odesolver( @(T, Y) m_BatchCultModel_DC(T, Y, hPR, xPR, y0scales, SysTopol), [max(t1) + (1e-6), tmax], Y0, odeoptions);
    
    % --- Indicate the switch has happened --------------------------------
    if y2(end,34) > 0 % iI is internal and so the switch happened ok
        switchInd = 1;
    else
        switchInd = 0;
    end
    
else % --- If no substrate then termination simualtion --------------------
    
    t2 = []; y2 = []; switchInd = 0;
    
end

%% ===== TEST QUALITY OF SIMULATION =======================================

%% ===== CALCULATE OUTPUTS TO RETURN ======================================
% --- Full time course ----------------------------------------------------
tsol = [t1; t2]; ysol = [y1; y2];

% --- Now unscale T and Y ---------------------------------------------
N0  = y0scales(1);
xS0 = y0scales(2);
xI0 = y0scales(3);
% Unscaled Y = scaledY * scaling vector
scalevector = ones(1, length(Y0));
% scalevector([1 2 3 4 5 6 27 34]) = [N0 xS0 xS0 xI0 xS0 xS0 xS0 xI0];
% scalevector([1 2 3 4]) = [N0 xS0 xS0 xI0];
Y = scalevector.*ysol;

% --- Merge time vector -----------------------------------------------
T = [t1; t2];

% --- Measuring production performance --------------------------------
if switchInd == 1
    
    % --- Volumetric productivity -------------------------------------
    %     =(total prod / total batch culture time)
    vProd  = Y(end,3)/T(end);
    
    % --- Product yield -----------------------------------------------
    %     = (total prod / total potential product)
    if SysTopol(1) == 1
        pYield = Y(end,3)/(Y(1,2));
    elseif SysTopol(1) == 2
        pYield = Y(end,3)/(hPR(2)*Y(1,2));
    end
    
    if pYield > 1
        ProdPerf = [NaN, NaN];
	flag = 0;
    else
    	ProdPerf = [vProd, pYield];
	flag = 1;
    end
        
else
    ProdPerf = [0, 0];
    flag = 0;
end

end

