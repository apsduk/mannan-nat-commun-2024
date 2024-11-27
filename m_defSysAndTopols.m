% ------------------------------------------------------------------------
% Make matrix of all topologies to consider for optimization
% Ahmad A. Mannan
% Edits Alexander P.S. Darlington
% ------------------------------------------------------------------------

function SysTopols = m_defSysAndTopols


% Defining the inducible control systems:
SysTopols = [];
regs = [-1,0,1]; % regulation modes (inhibition (-1), const ENGINEERED exp (0), activation (+1))
count = 0;
for pre = 1:2 % precursor 'is' (1) or 'ee' (2)
    for regT = [8, regs]  % TF reg of T, regulation is either host (8) or eng. regs [-1, 0 1]
        for regE = regs  % TF reg of E, regulation is eng. regs [-1, 0 1]
            for regEp = [-1,0]  % TF reg of Ep is eng. and either inhibited or unregulated
                for regTp =  [9, regs]  % TF reg of Tp is eng. regs [-1, 0 1] or absent 9
                    for regTF = regs  % TF reg of TF is eng. regs [-1, 0 1]
                        
                        % --- If Tp is absent then update xp
                        % exporter native 'X' alone (1) or supplemented with heterologous 'Tp' (2)
                        if regTp == 9
                            xp = 1;
                        else
                            xp = 2;
                        end
                                               
                        % --- If all genes are consitutive then skip all
                        if sum([regT regE regEp regTp] == [0 0 0 0]) == 4
                            % SKIP
                        else
                            
                            % --- Define possible topologies ------------------
                            count = count + 1;
                            SysTopols(count,:) = [pre,xp,regT,regE,regEp,regTp,regTF]; %#ok<AGROW>
                            
                        end
                        
                    end % regTF
                end % regTp
            end % regEp
        end % regE
    end % regT
end % pre

end
