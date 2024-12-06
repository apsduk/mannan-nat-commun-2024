% ------------------------------------------------------------------------
% TX regulation by inducible TF
% ------------------------------------------------------------------------

function r_g2m = m_TXreg(params,pTF)

% Define params:
eta    = params(1); % defines nature of TF regulation (0 = const exp; 1 = TF activated, -1 = TF inhibited)
w0     = params(2); % leaky expression as a fraction of maximum expression
g2m_mx = params(3); % maximum expression rate
K      = params(4); % affinity of TF for gene operator site

% Calculate transcription rate:
r_g2m = ( g2m_mx * (1 - eta^2) ) + (eta^2)*( (w0*g2m_mx) +  (g2m_mx/(1 +  K*pTF))*((K*pTF)^((1+eta)/2)) );
