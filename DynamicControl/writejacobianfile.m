
function [fname, J] = writejacobianfile(odesourcename, fname, T, Y, hPR, xPR, y0scales, SysTopol)

%% %%%%% MAKE SYMBOLIC VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sym_T = sym('T');
sym_Y = sym('Y', size(Y));
sym_hPR = sym('hPR', size(hPR));
sym_xPR = sym('xPR', size(xPR));
sym_y0scales = sym('y0scales', size(y0scales));
odefun = str2func(odesourcename);

%% %%%%% CALCULATE ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating symbolic derivative.');
dY_by_dt = odefun(sym_T, sym_Y, sym_hPR, sym_xPR, sym_y0scales, SysTopol);

%% %%%%% CALCULATE JACOBIAN SYMBOLICALLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating symbolic Jacobian.');
J = jacobian(dY_by_dt, sym_Y);

%% %%%%% WRITE JACOBIAN AS A FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Create jacobian name ------------------------------------------------
fname = [fname,'_autogen'];

disp(['Writing ',fname]);
[nrow, ncol] = size(J);

% --- Initilise new file --------------------------------------------------
fid = fopen([fname,'.m'],'w');
fprintf(fid, 'function J = jacode(T, Y, hPR, xPR, y0scales)');
fprintf(fid,'\n \n');

% --- Write species vectors -----------------------------------------------
fprintf(fid,'%% ===== Species ========== ========== ========== ==========\n');
for i = 1:numel(sym_Y)
    fprintf(fid, strcat(char(sym_Y(i)),' = Y(',num2str(i),'); \n'));
end
fprintf(fid,'\n');

% --- Write y0scales vectors -----------------------------------------------
fprintf(fid,'%% ===== Species ========== ========== ========== ==========\n');
for i = 1:numel(sym_y0scales)
    fprintf(fid, strcat(char(sym_y0scales(i)),' = y0scales(',num2str(i),'); \n'));
end
fprintf(fid,'\n');

% --- Write host parameters -----------------------------------------------
fprintf(fid,'%% ===== Host parameters ========== ========== ========== ==========\n');
for i = 1:numel(sym_hPR)
    fprintf(fid, strcat(char(sym_hPR(i)),' = hPR(',num2str(i),'); \n'));
end
fprintf(fid,'\n');

% --- Write circuit parameters --------------------------------------------
fprintf(fid,'%% ===== Circuit parameters ========== ========== ========== ==========\n');
for i = 1:numel(sym_xPR)
    fprintf(fid, strcat(char(sym_xPR(i)),' = xPR(',num2str(i),'); \n'));
end
fprintf(fid,'\n');

% --- Write Jacobian one element at a time --------------------------------
fprintf(fid,'%% ===== Jacobian ========== ========== ========== ==========\n');
fprintf(fid,['J = zeros(',num2str(nrow),',',num2str(ncol),'); \n']);
for j_1 = 1:length(J(1,:))
    for j_2 = 1:length(J(:,1))
        fprintf(fid, strcat('J(',num2str(j_1),',',num2str(j_2),') = ',char(J(j_1,j_2)),'; \n'));
    end
end
fprintf(fid,'\n');

% --- Make J sparse -------------------------------------------------------
% fprintf(fid, 'J = sparse(J);');

% --- Close and save output -----------------------------------------------
fprintf(fid, '\n end \n');
fclose(fid);
disp([fname,' file complete.']);

end

