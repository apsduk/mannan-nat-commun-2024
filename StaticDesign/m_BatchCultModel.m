% ------------------------------------------------------------------------
% BATCH CULTURE MODEL
% ------------------------------------------------------------------------

function [dY,lambda] = m_BatchCultModel(t,xvals,hPR,xPR,SysTopol)


% --- Define variables ---------------------------------------------------

X  = xvals(1);    % total biomass
xS = xvals(2);    % substrate in media
xP = xvals(3);    %#ok<*NASGU> % product in media
Ix = xvals(4);    % inducer in media
y = xvals(5:end); % all internal species


% --- Define substrate and inducer in media ------------------------------

% Define amount of substrate in media:
hPR_t = hPR;
if xS <= 0
    xS = 0;
end
% disp(xS)
hPR_t(1) = xS;

% Define amount of inducer in media:
xPR_t = xPR;
xPR_t(26) = Ix;


% --- Simulate dynamics of internal species ------------------------------

% Calculate the derivatives of all internal species:
[dy,lambda,~,~,rfluxes] = model_full(t,y,hPR_t,xPR_t,SysTopol);


% --- Simulate dynamics of external species and culture process ----------

% define reaction rates:
rT     = rfluxes(1); % substrate uptake rate (per cell)
r_uInd = rfluxes(6); % inducer uptake/secretion rate (per cell)
r_p    = rfluxes(5); % product secretion rate (per cell)

% total biomass:
dX  = (lambda * X);

% total substrate in the media:
if xS == 0
    dxS = 0;
else
    dxS = - (rT * X);
end

% total product in media:
dxP = (r_p * X);

% total inducer in media:
dIx = - (r_uInd * X);

% Vector of all derivative values:
dY = [dX; dxS; dxP; dIx; dy];