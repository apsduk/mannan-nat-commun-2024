% ------------------------------------------------------------------------
% Run model to steady state
% ------------------------------------------------------------------------

function [SSy,growth,fluxes,protmass,txrates,tlrates] = m_steadystate(hPR,xPR,y0,SysTopol)


% Run to solve, till within tolerance of steady state:
tol = 1e-5;
maxIters = 10000;
tspan = [0,10000];
i = 0;
while i <= maxIters
    % display iteration number
    i = i + 1;
    
    % solve ODEs:
    [t,y] = ode15s(@(t,y)model_full(t,y,hPR,xPR,SysTopol),tspan,y0);
    
    % check if steady state achieved:
    [dy,growth,protmass,~,fluxes,txrates,tlrates] = model_full(t,y(end,:),hPR,xPR,SysTopol);
    if max(abs(dy(1:19))) >= tol % ignoring possible accumulation of internal product (ip)
        y0 = y(end,:);
    else
        break
    end
end

if max(abs(dy(1:19))) >= tol % tell me if the system did not converge to steady state
    disp(['System did not converge to steady state! Max dy/dt = ', num2str(max(abs(dy(1:19))))])
end

% Save steady state of all intracellular species:
SSy = y(end,:);
