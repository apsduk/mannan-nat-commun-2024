% ------------------------------------------------------------------------
% Make matrix of all topologies to consider for optimization
% ------------------------------------------------------------------------

function SysTopols = m_defSysAndTopols


% Defining the inducible control systems:
SysTopols = [];
regs = [-1,0,1]; % regulation modes (inhibition (-1), const exp (0), activation (+1))
count = 0;
for p = 1:2 % precursor 'is' (1) or 'ee' (2)
    for xp = 1:2 % exporter native 'X' alone (1) or supplemented with heterologous 'Tp' (2)
        for regT = regs  % TF reg of T
            for regE = regs  % TF reg of E
                for regX = 0  % TF reg of X
                    for regEp = [-1,0]  % TF reg of Ep
                        for regTp = regs  % TF reg of Tp
                            for regTF = regs  % TF reg of TF
                                
                                % defining possible topologies:
                                if (p == 1 && regE ~= -1) || p == 2
                                    if xp == 1 
                                        if regTp == 0
                                            count = count + 1;
                                            SysTopols(count,:) = [p,xp,regT,regE,regX,regEp,regTp,regTF]; %#ok<AGROW>
                                        end
                                    elseif xp == 2
                                        count = count + 1;
                                        SysTopols(count,:) = [p,xp,regT,regE,regX,regEp,regTp,regTF]; %#ok<AGROW>
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
        end
    end
end
