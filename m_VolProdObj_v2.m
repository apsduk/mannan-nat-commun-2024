% ------------------------------------------------------------------------
% CALCULATING OBJECTIVES (I)  VOLUMETRIC PRODUCTIVITY
% ------------------------------------------------------------------------
% Alexander P.S. Darlington
% ------------------------------------------------------------------------
function [Objs, T, Y, flag, tsol, ysol] = m_VolProdObj(p, hPR0, xPR0, Y0scales, SysTopol, tmax, xidx, odesolver, odeoptions)

% --- Define parameters --------------------------------------------------

% params of host:
hPR = hPR0;

% params of circuit and synthesis pathway:
xPR = xPR0;

% update parameters in position xidx with the varying parameters p
xPR(xidx) = p(1:end-1);

% time of induction (t_ind):
tind  = p(end);

% --- Simulate batch culture ---------------------------------------------
idx = 0;
flag = 0;
while flag == 0
    idx = idx + 1;
    [Objs, T, Y, flag, tsol, ysol] = m_SimBatchCult_ProdPerf(hPR, xPR, Y0scales, SysTopol, tind, tmax, odesolver, odeoptions{idx});
    if idx == length(odeoptions)
        flag = 1;
    end
    % fprintf('odetolidx %d, flag %d, Obj %d \n', idx, flag, Objs(2));
end
Objs = (-1)*Objs(1);

end