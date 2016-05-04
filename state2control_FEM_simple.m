function control = state2control_FEM(state,i_a,grids,param)
% state = A, k, e
% coeff_me/mk are all coefficients
% param a container of parameters

% load parameters
bbeta = param.bbeta;
kkappa_S = param.kkappa_S;
kkappa_F = param.kkappa_F;
ttau = param.ttau;
h = param.h;
iiota = param.iiota;
aalpha = param.aalpha;
gyratio = param.gyratio;
ddelta = param.ddelta;

% load grids
Knodes = grids.Knodes;
Enodes = grids.Enodes;
EMKval = grids.EMKval;
EMEval = grids.EMEval;

A = state(1); k = state(2); e = state(3);

% find expectation terms
EMK = globaleval(k,e,Knodes,Enodes,squeeze(EMKval(i_a,:,:)));
EME = globaleval(k,e,Knodes,Enodes,squeeze(EMEval(i_a,:,:)));

% Find expectation terms tomorrow
Eme = EME; Emk = EMK;

% Back out other variables
c = (bbeta*Emk)^(-1);
q = kkappa_S/ttau*(1-c*bbeta*Eme-kkappa_F/ttau)^(-1);
if q >= 1
	q = 1.0;
	% warning('q>=1');
elseif q <=0
	q = 1e-7;
	% warning('q<=0');
end
ttheta = (q^(-iiota) - 1)^(1/iiota);
if ttheta ~= 0
	f = (1+ttheta^(-iiota))^(-1/iiota);
else
	f = 0;
end
U = 1-(kkappa_S/q+kkappa_F)*(f+(1-f)/ttau)-h;

y = A*k^aalpha;
v = ttheta*(e);
m = v*(e)/(v^iiota+(e)^iiota)^(1/iiota);

% Find mk and me;
mk = (1-ddelta+U*aalpha*y/k)/c;
me = U/c;

% Find next state
eplus = e + y - m;
I = y-c+e-eplus-kkappa_S*v-kkappa_F*m-e*h;
kplus = (1-ddelta)*k + I;

% pack stuff into control
control.kplus = kplus;
control.eplus = eplus;
control.mk = mk;
control.me = me;
control.q = q;

end