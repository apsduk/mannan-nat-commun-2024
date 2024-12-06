    % ------------------------------------------------------------------------
% MODEL OF HOST-AWARE PRODUCT SYNTHESIS - MOST GENERAL FORM
% ------------------------------------------------------------------------

function [dy,lambda,protmass,ribomassfrac,fluxes,TXrates,TLrates,summass] = model_full(t,y,hPR,xPR,SysTopol) %#ok<INUSL>


%% -- GET VARIABLES ------------------------------------------------------

y(y < 0) = 0;

% --- Define host variables ----------------------------------------------
% Note m = mRNA, c = complex, p = protein
is  = y( 1); ee  = y( 2);             % internal substrate (is) and energy (ee)
mT  = y( 3); cT  = y( 4); pT = y( 5); % transporter
mE  = y( 6); cE  = y( 7); pE = y( 8); % enzymes
mH  = y( 9); cH  = y(10); pH = y(11); % approx fixed amount of host proteins (q-fraction)
mR  = y(12); cR  = y(13); xR = y(14); rr = y(15); pR = y(16); % ribosomes
mX  = y(17); cX  = y(18); pX = y(19); % export reaction

% --- Define pathway variables -------------------------------------------
mEp = y(20); cEp = y(21); pEp = y(22);
ip  = y(23);                           % synthesis pathway enzyme
mTp = y(24); cTp = y(25); pTp = y(26); % product synthesis transporter
mTF = y(27); cTF = y(28); pTF = y(29); % inducible biosensor for synthetic inducer
I   = y(30); sC  = y(31);              % inducer and sequestered complex


%% -- GET PARAMETERS -----------------------------------------------------

% --- Define host system paramters ---------------------------------------
xS0     = hPR(1);  sS   = hPR(2);
vT      = hPR(3);  vE   = hPR(4);
KmT     = hPR(5);  KmE  = hPR(6);
wX      = hPR(7);  wH   = hPR(8);  wR = hPR(9); wr = hPR(10); % max transcription rates
oX      = hPR(11); oR   = hPR(12);
nX      = hPR(13); nR   = hPR(14);
bX      = hPR(15); uX   = hPR(16);
brho    = hPR(17); urho = hPR(18);
deg_m   = hPR(19);
kH      = hPR(20); hH   = hPR(21);
maxG    = hPR(22); kG   = hPR(23); M0 = hPR(24);
xphi    = hPR(25); vX   = hPR(26); KmX = hPR(27);

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
sTX_X   = xPR(11); % scaling maximum TX rate for native X
sTX_Ep  = xPR(12); % scaling maximum TX rate of Ep
sTX_Tp  = xPR(13); % scaling maximum TX rate of Tp
sTX_TF  = xPR(14); % scaling maximum TX rate of TF
K_T     = xPR(15); % biosensor affinity for native T
K_E     = xPR(16); % biosensor affinity for native E
K_X     = xPR(17); % biosensor afffinity for native exporter X
K_Ep    = xPR(18); % biosensor affinity for TX of Ep
K_Tp    = xPR(19); % biosensor affinity for TX of Tp
K_TF    = xPR(20); % biosensor affinity for TF PAR

% --- Define inducer parameters ------------------------------------------
kdiffI  = xPR(21); % rate of diffusion of inducer into and out of a single cell
VolCell = xPR(22); % volume of cell in L
VolCult = xPR(23); % working volume of culture in L
ksf     = xPR(24); % forward rate of TF sequestration by inducer I
ksr     = xPR(25); % reverse rate of TF sequestration by inducer I
Ix      = xPR(26); % inducer amount (molecules) added to working culture volume

% --- Define metabolic system --------------------------------------------
prodPre = SysTopol(1); % product precursor, 'is' (1) or 'ee' (2)
prodExp = SysTopol(2); % product exporter, 'nat' (1) or with 'het' (2), i.e. native alone or with heterologous exporter

% --- Define circuit topologies ------------------------------------------
ctT     = SysTopol(3); % indicator of TF control on TX of T
ctE     = SysTopol(4); % indicator of TF control on TX of E
ctX     = SysTopol(5); % indicator of TF control on TX of X
ctEp    = SysTopol(6); % indicator of TF control on TX of Ep
ctTp    = SysTopol(7); % indicator of TF control on TX of Tp
ctTF    = SysTopol(8); % indicator of TF autoregulation


%% -- DEFINE MODEL REACTION RATES ----------------------------------------

% --- Define transcription rates -----------------------------------------
% ee-dependent transcription rates with scaling
g2mH     = ((wH*ee)/(oX + ee))*(1/(1+(pH/kH)^hH));
g2mR     = ((wR*ee)/(oR + ee));
g2rr     = ((wr*ee)/(oR + ee));
g2mT_mx  = sTX_T  * ((wX*ee)/(oX + ee));
g2mE_mx  = sTX_E  * ((wX*ee)/(oX + ee));
g2mX_mx  = sTX_X  * ((xphi*wX*ee)/(oX + ee));
g2mEp_mx = sTX_Ep * ((wEp*ee)/(oX + ee));
g2mTp_mx = sTX_Tp * ((wTp*ee)/(oX + ee));
g2mTF_mx = sTX_TF * ((wTF*ee)/(oX + ee));

% include TF regulation of TX:
g2mT  = m_TXreg([ctT,  w0, g2mT_mx,  K_T],  pTF);
g2mE  = m_TXreg([ctE,  w0, g2mE_mx,  K_E],  pTF);
g2mX  = m_TXreg([ctX,  w0, g2mX_mx,  K_X],  pTF);
g2mEp = m_TXreg([ctEp, w0, g2mEp_mx, K_Ep], pTF);
g2mTp = m_TXreg([ctTp, w0, g2mTp_mx, K_Tp], pTF);
g2mTF = m_TXreg([ctTF, w0, g2mTF_mx, K_TF], pTF);

% --- Define translation rates -------------------------------------------
% Global translation rate (elongation rate):
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

% --- Growth rate --------------------------------------------------------
lambda = (1/M0)*gammaX*(cH + cR + cT + cE + cX + cEp + cTp + cTF);         % ===== UPDATE IF NEW GENES ADDED ====================

% --- Define Metabolic Reaction Rates ------------------------------------
rT = (pT*vT*xS0)/(KmT + xS0); % Substrate uptake rate
rE = (pE*vE*is)/(KmE + is);   % Host metabolism to biosynthetic precursor

% Production synthesis rate ...
if prodPre == 1     % ... consuming metabolite 'is'
    r_Psynth_is = (pEp*k_Ep*is)/(Km_Ep + is);
    r_Psynth_ee = 0;
elseif prodPre == 2 % ... consuming metabolite 'ee'
    r_Psynth_is = 0;
    r_Psynth_ee = (pEp*k_Ep*ee)/(Km_Ep + ee);
end

% Product secretion rate ...
r_X  = (pX*vX*ip)/(KmX + ip);
if prodExp == 1      % ... using native exporter
    r_Tp = 0;
elseif prodExp == 2  % ... supplementing with heterologous exporter
    r_Tp = (pTp*k_Tp*ip)/(Km_Tp + ip);
end
r_Exp = r_X + r_Tp;

% Inducer uptake rate:
r_uInd = kdiffI*(Ix*(VolCell/VolCult) - I);


%% -- DEFINE ODEs OF THE HOST --------------------------------------------

% Host metabolism:
dis = rT - rE - lambda*is - r_Psynth_is;                                   % ===== UPDATE IF NEW is CONSUMING REACTION ADDED ====================
dee = (sS*rE)  - lambda*ee - r_Psynth_ee ...                               % ===== UPDATE IF NEW ee CONSUMING OR NEW GENES REACTION ADDED ====================
     - nR*m2xR  - nX*m2pT   - nX*m2pE - nX*m2pX - nX*m2pH ...
     - nX*m2pEp - nX*m2pTp  - nX*m2pTF;

% Substrate transporter (T):
dmT = g2mT - (lambda + deg_m)*mT + m2pT - bX*pR*mT + uX*cT;
dcT = - lambda*cT + bX*pR*mT - uX*cT - m2pT;
dpT = m2pT - lambda*pT;

% Metabolic enzyme (E):
dmE = g2mE - (lambda + deg_m)*mE + m2pE - bX*pR*mE + uX*cE;
dcE = - lambda*cE + bX*pR*mE - uX*cE - m2pE;
dpE = m2pE - lambda*pE;

% House-keeping proteins (H):
dmH = g2mH - (lambda + deg_m)*mH + m2pH - bX*pR*mH + uX*cH;
dcH = - lambda*cH + bX*pR*mH - uX*cH - m2pH;
dpH = m2pH - lambda*pH;

% Native product exporter (X):
dmX = g2mX - (lambda + deg_m)*mX + m2pX - bX*pR*mX + uX*cX;
dcX = - lambda*cX + bX*pR*mX - uX*cX - m2pX;
dpX = m2pX - lambda*pX;

% Inactive ribosomes (xR):
dmR = g2mR - (lambda + deg_m)*mR + m2xR - bX*pR*mR + uX*cR;
dcR = - lambda*cR + bX*pR*mR - uX*cR - m2xR;
dxR = m2xR - lambda*xR - brho*xR*rr + urho*pR;

% Activated ribosome (in complex with ribosomal RNA), pR:
drr = g2rr - lambda*rr - brho*xR*rr + urho*pR;
dpR = brho*xR*rr - urho*pR - lambda*pR ...                                 % ===== UPDATE IF NEW GENES REACTION ADDED ====================
    + m2pT  - bX*pR*mT  + uX*cT ...
    + m2pE  - bX*pR*mE  + uX*cE ...
    + m2pH  - bX*pR*mH  + uX*cH ...
    + m2pX  - bX*pR*mX  + uX*cX ...
    + m2xR  - bX*pR*mR  + uX*cR ...
    + m2pEp - bX*pR*mEp + uX*cEp ...
    + m2pTp - bX*pR*mTp + uX*cTp ...
    + m2pTF - bX*pR*mTF + uX*cTF;


%% -- DEFINE ODEs OF THE SYNTHETIC PATHWAY -------------------------------

% Expression of production pathway enzyme (Ep):
dmEp = g2mEp - (lambda + deg_m)*mEp + m2pEp - bX*pR*mEp + uX*cEp;
dcEp = - lambda*cEp + bX*pR*mEp - uX*cEp - m2pEp;
dpEp = m2pEp - lambda*pEp;

% Metabolism - formation of intracellular product:
dip  = r_Psynth_is + r_Psynth_ee - r_Exp - lambda*ip;

% Expression of production pathway exporter enzyme (Tp):
dmTp = g2mTp - (lambda + deg_m)*mTp + m2pTp - bX*pR*mTp + uX*cTp;
dcTp = - lambda*cTp + bX*pR*mTp - uX*cTp - m2pTp;
dpTp = m2pTp - lambda*pTp;

% Expression of TF biosensor:
dmTF = g2mTF - (lambda + deg_m)*mTF + m2pTF - bX*pR*mTF + uX*cTF;
dcTF = - lambda*cTF + bX*pR*mTF - uX*cTF - m2pTF;
dpTF = m2pTF - lambda*pTF - (ksf*I^2*pTF) + (ksr*sC);
dI   = r_uInd - 2*((ksf*I^2*pTF) - (ksr*sC)) - lambda*I;
dCI  = (ksf*I^2*pTF) - (ksr*sC) - lambda*sC;


%% -- RETURN OUTPUTS -----------------------------------------------------

% Derivatives:
dy = [dis dee ...
      dmT dcT dpT ...
      dmE dcE dpE ...
      dmH dcH dpH ...
      dmR dcR dxR drr dpR ...
      dmX dcX dpX ...
      dmEp dcEp dpEp ...
      dip ...
      dmTp dcTp dpTp ...
      dmTF dcTF dpTF ...
      dI dCI]';
% TX and TL rates:
TXrates = [g2mT g2mE g2mH g2mR g2rr g2mX g2mEp g2mTp g2mTF];
TLrates = [m2pT m2pE m2pH m2xR m2pX m2pEp m2pTp m2pTF];

% Protein masses and ribosomal mass fraction:
protmass = [nX*pT, nX*pE, nX*pX, nX*pH, nR*(xR + pR + cT + cE + cX + cH + cR + cEp + cTp + cTF), nX*pEp, nX*pTp, nX*(pTF + sC)];
summass = sum(protmass);
ribomassfrac = (nR/M0)*(xR + pR + cT + cE + cX + cH + cR + cEp + cTp + cTF);

% Fluxes:
fluxes = [rT rE r_Psynth_is r_Psynth_ee r_Exp r_X r_Tp r_uInd];

