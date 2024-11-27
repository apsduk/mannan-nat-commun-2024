% ------------------------------------------------------------------------
% Defining model parameters
% Edits Alexander P.S. Darlington
% ------------------------------------------------------------------------

function [h_params,x_params] = model_params(sS0, vX0, KmX0, VolCult0)


% --- Parameters of host components --------------------------------------
xS0     = 1e4;          % 1 extracellular substrate concentration
sS      = sS0;          % 2 nutrient efficiency (i.e. stoichiometry of substrate to precursors conversion)
vT      = 728;          % 3 kcat of transport rxn
vE      = 5800;         % 4 kcat of metabolic rxn (with enzyme)
KmT     = 1e3;          % 5 Km of transport rxn
KmE     = 1e3;          % 6 Km of metabolic rxn
wX      = 4.14;         % 7 max transcription rate of metabolic proteins (transporter and enzymes)
wH      = 948.93;       % 8 max transcription rate of house-keeping proteins
wR      = 930;          % 9 max transcription rate of ribosomal protein
wr      = 3170;         % 10 max transcription rate of ribosomal RNA
oX      = 4.38;         % 11 effective Km of transcription of mRNA of metabolic proteins (T,H,E)
oR      = 426.87;       % 12 effective Km of transcription of ribosomal mRNA and RNA
nX      = 300;          % 13 average number of codons on mRNA of metabolic proteins (T,H,E), for translation
nR      = 7459;         % 14 average number of codons on mRNA of ribosomal mRNA, for translation
bX      = 1;            % 15 forward rate of active ribosomes to form complex with mRNA
uX      = 1;            % 16 rate of unbinding ribosome from mRNA
brho    = 1;            % 17 rate of binding ribosomal RNA to ribosome protein part
urho    = 1;            % 18 rate of unbinding of ribosomal RNA to ribosome protein part
deg_m   = 0.1;          % 19 average mRNA degradation rate
kH      = 152219 * 0.8; % 20 inverse threshold of conc of house-keeping proteins to inhibit their own expression (phenomenological - to maintain fixed levels for diff conditions)
hH      = 4 * 2;        % 21 Hill coefficient of house-keeping proteins to inhibit their own transcriptional expression
maxG    = 1260;         % 22 maximal translation rate
kG      = 7;            % 23 threshold conc of precursor for half-max translation rate
M0      = 1e8;          % 24 total cell mass

% Parameters of promiscuous host exporter:
xphi    = 0.05;  % 25 scaling factor for max transcription rate of host co-opted transporter X, TX = xphi*wX
vX      = vX0;    % 26 kcat of ip export rxn
KmX     = KmX0;   % 27 Km of ip export rxn


% --- Parameters of exogenous components ---------------------------------

% Parameters of TX expression and kinetics of heterologous enzymes:
w0      = 1e-4;            % 1 leakiness of the artifical promoters as a fraction of their maximum transcription rate (%)
wEp     = 20;              % 2 max transcription rate of Ep
k_Ep    = vE; Km_Ep = KmE; % 5,6 kcat and Km for kinetics of Ep
wTp     = 20;              % 3 max transcription rate of Tp
k_Tp    = vT; Km_Tp = KmT; % 1e5*KmT; % 7,8 kcat and Km for kinetics of Tp

% Parameters of inducible control system:
wTF     = 20;              % 4 max transcription rate of TF

% Parameters of the control system
sTX_T   = 1;               %  9 scaling maximum TX rate of native T
sTX_E   = 1;               % 10 scaling maximum TX rate of native E
sTX_Ep  = 0;               % 11 scaling maximum TX rate of Ep
sTX_Tp  = 0;               % 12 scaling maximum TX rate of Tp
sTX_TF  = 0;               % 13 scaling maximum TX rate of TF

K_T     = 1/(1000/300);    % 14 (1/molecules) - biosensor affinity for native T
K_E     = 1/(1000/300);    % 15 (1/molecules) - biosensor affinity for native E
K_Ep    = 1/(1000/300);    % 16 (1/molecules) assuming 1nM = 1 molecule - (1000/(300 1/uM)) - based on affinity of 300uM from Mannan & Bates (2021) - biosensor affinity to bind to inhibit TX of Ep
K_Tp    = 1/(1000/300);    % 17 (1/molecules) - biosensor affinity for TX of Tp
K_TF    = 1/(1000/300);    % 18 (1/molecules) - biosensor affinity for TF PAR

% Parameters for inducer uptake and induction:
kdiffI  = 3600/60;         % 19 1/min - rate of diffusion of inducer into and out of a single cell
VolCell = 1e-15;           % 20 volume of cell in L - based on 1uM wide by 2uM sized E. coli cell (Neidhardt (1990))
VolCult = VolCult0;        % 21 working volume of culture in L - based on a 3L benchtop vessel
ksf     = 1057.7/60;       % 22 1/(molecules^2.min) - forward rate of TF sequestration by inducer I - based on Mannan & Bates (2021)
ksr     = 1292.1/60;       % 23 1/min - reverse rate of TF sequestration by inducer I, based on Mannan & Bates (2021)

% Inducer amount added to working culture volume:
xI0     = 0;               % 24 (# of molecules)

% --- Return vectors -----------------------------------------------------

% Define vector of ...
% ... params of endogenous components:
% h_params = [xS0 sS vT vE KmT KmE wX wH wR wr oX oR nX nR bX uX brho urho deg_m kH hH maxG kG M0 xphi vX KmX]';
h_params   = [xS0 sS vT vE KmT KmE wX wH wR wr oX oR nX nR bX uX brho urho deg_m kH hH maxG kG M0 xphi vX KmX]';

% ... params of exogenous components:
% x_params = [w0 wEp wTp wTF k_Ep Km_Ep k_Tp Km_Tp sTX_T sTX_E sTX_Ep sTX_Tp sTX_TF K_T K_E K_Ep K_Tp K_TF kdiffI VolCell VolCult ksf ksr xI]';
x_params   = [w0 wEp wTp wTF k_Ep Km_Ep k_Tp Km_Tp sTX_T sTX_E sTX_Ep sTX_Tp sTX_TF K_T K_E K_Ep K_Tp K_TF kdiffI VolCell VolCult ksf ksr xI0]';

end
