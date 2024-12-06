% ------------------------------------------------------------------------
% GROWTH AND PRODUCTION IN BATCH CULTURE
% ------------------------------------------------------------------------
% Ahmad A. Mannan
% Edits Alexander P.S. Darlington
% ------------------------------------------------------------------------

function [dY, lambda, TXrates, TLrates, gammaX, protmass, summass, ribomassfrac, fluxes] = m_BatchCultModel_DC(T, Y, hPR, xPR, y0scales, SysTopol)

%% ===== VARIABLES ========================================================

% --- Scaling factors -----------------------------------------------------
N0  = y0scales(1);
xS0 = y0scales(2);
xI0 = y0scales(3);

% --- Define external variables -------------------------------------------
N  = Y(1); % total biomass
xS = Y(2); % substrate in media
xP = Y(3); % product in media
xI = Y(4); % inducer conc in media

% --- Define host variables -----------------------------------------------
% Note m = mRNA, c = complex, p = protein
iS  = Y( 5); ee  = Y( 6);             % internal substrate (is) and energY (ee)
mT  = Y( 7); cT  = Y( 8); pT = Y( 9); % transporter
mE  = Y(10); cE  = Y(11); pE = Y(12); % enzymes
mH  = Y(13); cH  = Y(14); pH = Y(15); % approx fixed amount of host proteins (q-fraction)
mX  = Y(16); cX  = Y(17); pX = Y(18); % export reaction
mR  = Y(19); cR  = Y(20); xR = Y(21); rr = Y(22); pR = Y(23); % ribosomes

% --- Define pathway variables --------------------------------------------
mEp = Y(24); cEp = Y(25); pEp = Y(26); iP  = Y(27); % synthesis pathway enzyme
mTp = Y(28); cTp = Y(29); pTp = Y(30);              % product synthesis transporter
mTF = Y(31); cTF = Y(32); pTF = Y(33);              % inducible biosensor for synthetic inducer
iI  = Y(34); sC  = Y(35);                           % inducer and sequestered complex

%% ===== PARAMETERS =======================================================

% --- Define host sYstem paramters ---------------------------------------
% xS0     = hPR(1);
sS    = hPR( 2);
vT    = hPR( 3); vE   = hPR( 4);
KmT   = hPR( 5); KmE  = hPR( 6);
wX    = hPR( 7); wH   = hPR( 8); wR = hPR(9); wr = hPR(10);
oX    = hPR(11); oR   = hPR(12);
nX    = hPR(13); nR   = hPR(14);
bX    = hPR(15); uX   = hPR(16);
brho  = hPR(17); urho = hPR(18);
deg_m = hPR(19);
kH    = hPR(20); hH   = hPR(21);
maxG  = hPR(22); kG   = hPR(23); M0 = hPR(24);
xphi  = hPR(25); vX   = hPR(26); KmX = hPR(27);

% --- Define pathway parameters ------------------------------------------
w0      = xPR( 1);  % leakiness of synthetic regulated promoters
wEp     = xPR( 2);  % maximum TX rate of Ep
wTp     = xPR( 3);  % maximum TX rate of Tp
wTF     = xPR( 4);  % maximum TX rate of TF biosensor
k_Ep    = xPR( 5); Km_Ep = xPR( 6); % kcat and KM for Ep
k_Tp    = xPR( 7); Km_Tp = xPR( 8); % kcat and KM for Tp

% --- Define pathway regulation parameters -------------------------------
sTX_T   = xPR( 9); % scaling maximum TX rate of native T
sTX_E   = xPR(10); % scaling maximum TX rate of native E
sTX_Ep  = xPR(11); % scaling maximum TX rate of Ep
sTX_Tp  = xPR(12); % scaling maximum TX rate of Tp
sTX_TF  = xPR(13); % scaling maximum TX rate of TF
K_T     = xPR(14); % biosensor affinitY for native T
K_E     = xPR(15); % biosensor affinitY for native E
K_Ep    = xPR(16); % biosensor affinitY for TX of Ep
K_Tp    = xPR(17); % biosensor affinitY for TX of Tp
K_TF    = xPR(18); % biosensor affinitY for TF PAR

% --- Define inducer parameters ------------------------------------------
kdiffI  = xPR(19); % rate of diffusion of inducer into and out of a single cell
VolCell = xPR(20); % volume of cell in L
VolCult = xPR(21); % working volume of culture in L
ksf     = xPR(22); % forward rate of TF sequestration by inducer I
ksr     = xPR(23); % reverse rate of TF sequestration by inducer I
% xI0     = xPR(24); % inducer amount (molecules) added to working culture volume

% --- Define metabolic sYstem --------------------------------------------
prodPre = SysTopol(1); % product precursor, 'is' (1) or 'ee' (2)
prodExp = SysTopol(2); % product exporter, 'nat' (1) or with 'het' (2), i.e. native alone or with heterologous exporter

% --- Define circuit topologies ------------------------------------------
ctT     = SysTopol(3); % TF control on TX of T
ctE     = SysTopol(4); % TF control on TX of E
ctEp    = SysTopol(5); % TF control on TX of Ep
ctTp    = SysTopol(6); % TF control on TX of Tp
ctTF    = SysTopol(7); % TF autoregulation

%% ===== CALCULATE RATES ==================================================

% --- Define transcription rates ------------------------------------------
% ee-dependent transcription rates with scaling
g2mH     = ((wH*ee)/(oX + ee))*(1/(1+(pH/kH)^hH));
g2mR     = ((wR*ee)/(oR + ee));
g2rr     = ((wr*ee)/(oR + ee));
g2mX     = xphi*(wX*ee)/(oX + ee);
g2mT_mx  = sTX_T  * ((wX*ee)/(oX + ee));
g2mE_mx  = sTX_E  * ((wX*ee)/(oX + ee));
g2mEp_mx = sTX_Ep * ((wEp*ee)/(oX + ee));
g2mTp_mx = sTX_Tp * ((wTp*ee)/(oX + ee));
g2mTF_mx = sTX_TF * ((wTF*ee)/(oX + ee));

% include TF regulation of TX:
g2mT  = m_TXreg([ctT,  w0, g2mT_mx,  K_T],  pTF);
g2mE  = m_TXreg([ctE,  w0, g2mE_mx,  K_E],  pTF);
g2mEp = m_TXreg([ctEp, w0, g2mEp_mx, K_Ep], pTF);
g2mTp = m_TXreg([ctTp, w0, g2mTp_mx, K_Tp], pTF);
g2mTF = m_TXreg([ctTF, w0, g2mTF_mx, K_TF], pTF);

% --- Define translation rates --------------------------------------------
% global translation rate (elongation rate):
gammaX = (maxG*ee)/(kG + ee);

% protein per translation complex per min)
m2pT  = (gammaX/nX)*cT;
m2pE  = (gammaX/nX)*cE;
m2pH  = (gammaX/nX)*cH;
m2xR  = (gammaX/nR)*cR;
m2pX  = (gammaX/nX)*cX;
m2pTp = (gammaX/nX)*cTp;
m2pEp = (gammaX/nX)*cEp;
m2pTF = (gammaX/nX)*cTF;

% --- Growth rate ---------------------------------------------------------
lambda = (1/M0)*gammaX*(cH + cR + cT + cE + cX + cEp + cTp + cTF);         % ===== UPDATE IF NEW GENES ADDED ====================

% --- Define Metabolic Reaction Rates -------------------------------------
r_U = (pT*vT*xS)/(KmT + xS); % Substrate uptake rate
r_E = (pE*vE*iS)/(KmE + iS); % Host metabolism to biosynthetic precursor

% --- Production synthesis ------------------------------------------------
if prodPre == 1     % ... consuming metabolite 'is'
    r_Psynth_is = (pEp*k_Ep*iS)/(Km_Ep + iS);
    r_Psynth_ee = 0;
elseif prodPre == 2 % ... consuming metabolite 'ee'
    r_Psynth_is = 0;
    r_Psynth_ee = (pEp*k_Ep*ee)/(Km_Ep + ee);
end

% --- Product secretion rate ----------------------------------------------
if prodExp == -1
    r_Exp = r_Psynth_is + r_Psynth_ee;
elseif prodExp == 1      % ... using native exporter
    r_Exp = (pX*vX*iP)/(KmX + iP);
elseif prodExp == 2  % ... supplementing with heterologous exporter
    r_Exp = ((pX*vX*iP)/(KmX + iP)) + ((pTp*k_Tp*iP)/(Km_Tp + iP));
end

% --- Inducer uptake rate -------------------------------------------------
r_uInd = kdiffI*(xI*(VolCell/VolCult) - iI);

% --- Ensure all rates are > 0 --------------------------------------------
% if g2mT < 0; g2mT = 0; end
% if g2mE < 0; g2mE = 0; end
% if g2mEp < 0; g2mEp = 0; end
% if g2mTp < 0; g2mTp = 0; end
% if g2mTF  < 0; g2mTF = 0; end
% if gammaX < 0; gammaX = 0; end
% if m2pT < 0; m2pT = 0; end
% if m2pE < 0; m2pE = 0; end
% if m2pH < 0; m2pH = 0; end
% if m2xR < 0; m2xR = 0; end
% if m2pX < 0; m2pX = 0; end
% if m2pTp < 0; m2pTp = 0; end
% if m2pEp < 0; m2pEp = 0; end
% if m2pTF < 0; m2pTF = 0; end
% if lambda < 0; lambda = 0; end
% if r_U < 0; r_U = 0; end
% if r_E < 0; r_E = 0; end
% if r_Psynth_is < 0; r_Psynth_is = 0; end
% if r_Psynth_ee < 0; r_Psynth_ee = 0; end
% if r_Exp < 0; r_Exp = 0; end
% if r_uInd < 0; r_uInd = 0; end

%% ===== ENVIRONMENTAL ODEs ===============================================

% --- total biomass -------------------------------------------------------
dN = (lambda*N);

% --- total substrate in the media ----------------------------------------
% if xS < 0
%     dxS = 0;
% else
    dxS = - (r_U*N);
% end

% --- total product in media ----------------------------------------------
dxP = (r_Exp*N);

% --- total inducer in media ----------------------------------------------
dxI = - (r_uInd*N);

%% ===== HOST ODEs ========================================================
% --- host metabolism -----------------------------------------------------
diS = r_U - r_E - lambda*iS - r_Psynth_is;                                 % ===== UPDATE IF NEW is CONSUMING REACTION ADDED ====================
dee = (sS*r_E)  - lambda*ee - r_Psynth_ee ...                              % ===== UPDATE IF NEW ee CONSUMING OR NEW GENES REACTION ADDED ====================
    - nR*m2xR  - nX*m2pT  - nX*m2pE - nX*m2pX - nX*m2pH ...
    - nX*m2pEp - nX*m2pTp - nX*m2pTF;

% --- substrate transporter (T) -------------------------------------------
dmT = g2mT - (lambda + deg_m)*mT + m2pT - bX*pR*mT + uX*cT;
dcT = - lambda*cT + bX*pR*mT - uX*cT - m2pT;
dpT = m2pT - lambda*pT;

% --- metabolic enzyme (E) ------------------------------------------------
dmE = g2mE - (lambda + deg_m)*mE + m2pE - bX*pR*mE + uX*cE;
dcE = - lambda*cE + bX*pR*mE - uX*cE - m2pE;
dpE = m2pE - lambda*pE;

% --- house-keeping proteins (H) ------------------------------------------
dmH = g2mH - (lambda + deg_m)*mH + m2pH - bX*pR*mH + uX*cH;
dcH = - lambda*cH + bX*pR*mH - uX*cH - m2pH;
dpH = m2pH - lambda*pH;

% --- native product exporter (X) -----------------------------------------
dmX = g2mX - (lambda + deg_m)*mX + m2pX - bX*pR*mX + uX*cX;
dcX = - lambda*cX + bX*pR*mX - uX*cX - m2pX;
dpX = m2pX - lambda*pX;

% --- inactive ribosomes (xR) and rRNA (rr) -------------------------------
dmR = g2mR - (lambda + deg_m)*mR + m2xR - bX*pR*mR + uX*cR;
dcR = - lambda*cR + bX*pR*mR - uX*cR - m2xR;
dxR = m2xR - lambda*xR - brho*xR*rr + urho*pR;
drr = g2rr - lambda*rr - brho*xR*rr + urho*pR;

% -- activated ribosome (in complex with ribosomal RNA), pR ---------------
dpR = brho*xR*rr - urho*pR - lambda*pR ...                                 % ===== UPDATE IF NEW GENES REACTION ADDED ====================
    + m2pT  - bX*pR*mT  + uX*cT ...
    + m2pE  - bX*pR*mE  + uX*cE ...
    + m2pH  - bX*pR*mH  + uX*cH ...
    + m2pX  - bX*pR*mX  + uX*cX ...
    + m2xR  - bX*pR*mR  + uX*cR ...
    + m2pEp - bX*pR*mEp + uX*cEp ...
    + m2pTp - bX*pR*mTp + uX*cTp ...
    + m2pTF - bX*pR*mTF + uX*cTF;

%% ===== pathway AND CIRCUIT ODEs =========================================
% --- pathway enzyme (Ep) -------------------------------------------------
dmEp = g2mEp - (lambda + deg_m)*mEp + m2pEp - bX*pR*mEp + uX*cEp;
dcEp = - lambda*cEp + bX*pR*mEp - uX*cEp - m2pEp;
dpEp = m2pEp - lambda*pEp;

% --- metabolism - formation of intracellular product ---------------------
diP  = r_Psynth_is + r_Psynth_ee - r_Exp - lambda*iP;

% --- pathway exporter enzyme (Tp) ----------------------------------------
dmTp = g2mTp - (lambda + deg_m)*mTp + m2pTp - bX*pR*mTp + uX*cTp;
dcTp = - lambda*cTp + bX*pR*mTp - uX*cTp - m2pTp;
dpTp = m2pTp - lambda*pTp;

% --- TF biosensor --------------------------------------------------------
dmTF = g2mTF - (lambda + deg_m)*mTF + m2pTF - bX*pR*mTF + uX*cTF;
dcTF = - lambda*cTF + bX*pR*mTF - uX*cTF - m2pTF;
dpTF = m2pTF - lambda*pTF - ksf*(iI^2)*pTF + ksr*sC;
diI   = r_uInd - 2*ksf*(iI^2)*pTF + 2*ksr*sC - lambda*iI;
dsC  = ksf*(iI^2)*pTF - ksr*sC - lambda*sC;

%% ===== RETURN OUTPUTS ===================================================

% --- Update derivatives --------------------------------------------------
ddN  = dN; % /N0;
ddxS = dxS; % /xS0;
ddxP = dxP; % /xS0;
ddxI = dxI; % /xI0;
ddiS = diS; % /xS0;
ddee = dee; % /xS0;
ddiP = diP;
ddiI = diI; % /xI0;

% --- derivatives ---------------------------------------------------------
dY = [ddN; ddxS; ddxP; ddxI; ...
    ddiS; ddee; ...
    dmT; dcT; dpT; ...
    dmE; dcE; dpE; ...
    dmH; dcH; dpH; ...
    dmX; dcX; dpX; ...
    dmR; dcR; dxR; drr; dpR; ....
    dmEp; dcEp; dpEp; ...
    ddiP; ...
    dmTp; dcTp; dpTp;
    dmTF; dcTF; dpTF; ddiI; dsC];

% --- TX and TL rates -----------------------------------------------------
TXrates = [g2mT g2mE g2mH g2mR g2rr g2mX g2mEp g2mTp g2mTF];
TLrates = [m2pT m2pE m2pH m2xR m2pX m2pEp m2pTp m2pTF];

% --- protein masses and ribosomal mass fraction --------------------------
protmass = [nX*pT, nX*pE, nX*pX, nX*pH, nR*(xR + pR + cT + cE + cH + cR + cX + cEp + cTp + cTF), nX*pEp, nX*pTp, nX*(pTF + sC)];
summass = sum(protmass);
ribomassfrac = (nR/M0)*(xR + pR + cT + cE + cH + cR + cX + cEp + cTp + cTF);

% ---- fluxes -------------------------------------------------------------
fluxes = [r_U, r_E, r_Psynth_is, r_Psynth_ee, r_Exp, r_uInd];

end