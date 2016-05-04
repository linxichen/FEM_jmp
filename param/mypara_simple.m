% copied from dynare
aalpha  = 1/3;
bbeta   = 0.997;
ddelta  = 0.025;
uupsilon= 1.0;
Frisch  = 1.0;

% goods market frictions (calibration 1, qbar = 0.11)
options = optimoptions('fsolve','TolFun',1e-9);
kkappa_S = 0.01; % play with this
kkappa_F = 0.001; % play with this
ttau = 0.4; % share of surplus to retailer
iiota = 1.27;
h = 0.005;
ff = @(t) (1+t^(-iiota))^(-1/iiota);
qq = @(t) (1+t^(iiota))^(-1/iiota);
UU = @(t) (ff(t)*(1-kkappa_F-kkappa_S/qq(t))-h )/(1-bbeta*(1-ff(t)));
objfunc = @(t)  ttau*(1-bbeta*UU(t))-kkappa_S/qq(t)-kkappa_F;
tthetabar = fsolve(objfunc,100,options);
qbar = qq(tthetabar)
fbar = ff(tthetabar)

%% plot
figure
fplot(ff,[0 100])
hold on
fplot(qq,[0 100])

%%
% capital adjustment
adjcost = 0.0; 

% targets
cy      = 0.59;
nbar    = 0.3;
% fbar    = 0.4; % to match around 2.4 ISratio for final sale
% qbar    = 0.90;
gyratio = 0.15;
holdingcost = 0.04;

% Exo process
rrho_z      = 0.95;
rrho_zk     = 0.5;
rrho_g      = 0.95;
ssigma_z    = 0.008;
ssigma_zk   = 0.01;
ssigma_g    = 0.01;

% find SS and rest of parameters
%% Matching related
ey        = 1/fbar;
Pmbar = 1 - kkappa_S/qbar - kkappa_F;
Ubar = (1-kkappa_F/ttau-kkappa_S/ttau/qbar)/bbeta;

% Great ratios
kbar = ((1/bbeta-1+ddelta)/Ubar/aalpha)^(1/(aalpha-1));
rbar = Ubar*aalpha*kbar^(aalpha-1);
ybar      = kbar^aalpha;
ebar      = ey*ybar;
vbar      = ebar*tthetabar;
cbar      = (1-kkappa_F)*ybar - kkappa_S*vbar - ddelta*kbar - ebar*h;
mbar      = ybar;
zbar      = 1;
zkbar     = 1;
mmubar    = 1/cbar;

% pack into container
param.bbeta = bbeta;
param.kkappa_S = kkappa_S;
param.kkappa_F = kkappa_F;
param.ttau = ttau;
param.h = h;
param.iiota = iiota;
param.aalpha = aalpha;
param.gyratio = gyratio;
param.ddelta = ddelta;

% save 'param.mat';


