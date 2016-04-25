% copied from dynare
aalpha  = 0.33;
bbeta   = 0.997;
ddelta  = 0.025;
uupsilon= 1.0;
Frisch  = 1.0;

% goods market frictions (calibration 1, qbar = 0.11)
kkappa_S = 0.00000001; % play with this
kkappa_F = 0.16; % play with this
ttau = 0.2; % share of surplus to retailer

% capital adjustment
adjcost = 0.1; 

% targets
cy      = 0.59;
nbar    = 0.3;
fbar    = 0.4; % to match around 2.4 ISratio for final sale
qbar    = 0.90;
gyratio = 0.15;
holdingcost = 0.04;

% Exo process
rrho_z      = 0.95;
rrho_zk     = 0.5;
rrho_g      = 0.95;
ssigma_z    = 0.018;
ssigma_zk   = 0.01;
ssigma_g    = 0.01;

% find SS and rest of parameters
%% Matching related
ey        = 1/fbar - 1;
Pmbar = 1 - kkappa_S/qbar - kkappa_F;
term1 = -1 + (1-Pmbar)/ttau + bbeta*fbar*Pmbar/(1-bbeta*(1-fbar));
term2 = 1- bbeta*(1-fbar);
h = term1*term2;
tthetabar = fbar/qbar;
findela   = @(i) log(1+tthetabar^(-i)) + i*log(fbar);
options = optimoptions('fsolve','TolFun',1e-10);
iiota     = fsolve(findela,1.5,options);
Ubar      = (fbar*(1+h-kkappa_S/qbar-kkappa_F) - h) /(1-bbeta*(1-fbar));

% Great ratios
rbar      = 1/bbeta + ddelta - 1;
kovern    = (rbar/Ubar/aalpha)^(1/(aalpha-1));
ybar      = kovern^aalpha*nbar;
gbar      = gyratio*ybar;
ebar      = ey*ybar;
vovern    = tthetabar*(ebar+ybar)/nbar;
covern    = (1-kkappa_F)*ybar/nbar - kkappa_S*vovern - ddelta*kovern - ebar*h/nbar - gbar/nbar;

cbar      = covern*nbar;
kbar      = kovern*nbar;
ibar      = ddelta*kbar;
mbar      = ybar;
zbar      = 1;
zkbar     = 1;
mmubar    = 1/cbar;
Pmbar     = 1-(kkappa_S/qbar+kkappa_F);
gbar      = gyratio*ybar;
masterwedgebar = (1-aalpha)*ybar*(1-nbar)/nbar/cbar;


% Pick ppsi to satisfy cy = 0.6 and n = 0.2
wbar      = Ubar*(1-aalpha)*kovern^(aalpha);
pphi      = wbar*(1-nbar)/cbar;

% Find vacancy level from RC
vbar      = ((1-kkappa_F)*ybar - cbar - ddelta*kbar - ebar*h -gbar)/kkappa_S;

% Getting the rest
MCbar     = wbar^(1-aalpha)*rbar^(aalpha);


save param.mat;


