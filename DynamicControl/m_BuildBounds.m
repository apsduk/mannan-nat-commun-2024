% ------------------------------------------------------------------------
% BUILD BOUNDS BASED ON TOPOLOGY
% Alexander P.S.Darlington
% ------------------------------------------------------------------------

function [lb, ub, xidx, xnames, topodescr] = m_BuildBounds(SysTopol, sTX_Z, K_Z, t_ind)

% --- Indicies of key parameters in xPR -----------------------------------
sTX_idx(1) =  9; K_idx(1) = 14; parname{1} = 'T';
sTX_idx(2) = 10; K_idx(2) = 15; parname{2} = 'E';
sTX_idx(3) = 11; K_idx(3) = 16; parname{3} = 'Ep';
sTX_idx(4) = 12; K_idx(4) = 17; parname{4} = 'Tp';
sTX_idx(5) = 13; K_idx(5) = 18; parname{5} = 'TF';

% --- Initiate the topodescr ----------------------------------------------
if SysTopol(1) == 1
    topodescr = 'Substrate iS |';
elseif SysTopol(1) == 2
    topodescr = 'Substrate ee |';
end
if SysTopol(2) == 1
    topodescr = [topodescr,' export X |'];
elseif SysTopol(2) == 2
    topodescr = [topodescr,' export X + Tp |'];
end

% --- Build the bounds for the optimisation -------------------------------
fid = 0;
bnd = [];
xidx = [];
xnames = [];
fprintf('Building bounds for ... [ ');
for k = 1:length(SysTopol)-2
    
    % --- Print output ----------------------------------------------------
    if k == length(SysTopol)-2
        fprintf(num2str(SysTopol(k+2)));
    else
        fprintf([num2str(SysTopol(k+2)),', ']);
    end
    
    % --- Update the bounds and index for this gene TX rate ---------------
    if ismember(SysTopol(k+2), [-1 0 +1])
        
        % update parameters
        fid = fid+1;
        
        % update the bounds of the expression
        bnd(fid,:) = sTX_Z;
        xidx(fid) = sTX_idx(k);
        xnames{fid} = ['sTX_',parname{k}];
        
        % update description
        topodescr = [topodescr,' ',parname{k},'('];
        
        % --- Add regulation --------------------------------------------------
        if ismember(SysTopol(k+2),[ -1, +1])
            
            % update paraemter
            fid = fid + 1;
            
            % update the bounds and index for this gene regulation
            bnd(fid,:) = K_Z;
            xidx(fid) = K_idx(k);
            xnames{fid} = ['K_',parname{k}];
            
            % update description
            topodescr = [topodescr,num2str(SysTopol(k+2))];
            
        else
            
            % update description
            topodescr = [topodescr,'const.'];
            
        end
        
        % update description
        topodescr = [topodescr,')'];
        
    end
    
end

fprintf(']\n');

% --- Update bnds with t_ind ----------------------------------------------
xnames{end+1} = 't_ind';
bnd = [bnd; t_ind];
lb = bnd(:,1);
ub = bnd(:,2);

end