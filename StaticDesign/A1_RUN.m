% ------------------------------------------------------------------------
% RUN MULTI-OBJECTIVE OPTIMIZATION FOR STATIC DESIGNS
% ------------------------------------------------------------------------

function A1_RUN


%%
% [I] Multi-obj opt for static designs of ...
%     - host metabolite "is" and 
%     - unbound export from exporter X (no additional transporter Tp)

is = 1;
xporter = 2;
TuneEpAndE = 1;
fig = 1;
A1_MultiobjOpt(is,xporter,TuneEpAndE,fig)
