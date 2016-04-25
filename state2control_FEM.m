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
pphi = param.pphi;
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
q = kkappa_S/ttau*(1+h-c*bbeta*Eme-kkappa_F/ttau)^(-1);
if q >= 1
	q = 0.5;
	% warning('q>=1');
elseif q <=0
	q = 0.5;
	% warning('q<=0');
end
ttheta = (q^(-iiota) - 1)^(1/iiota);
f = (1+ttheta^(-iiota))^(1/iiota);
U = 1-(kkappa_S/q+kkappa_F)*(f+(1-f)/ttau);
GG = pphi*c/U/(1-aalpha)/A/(k^aalpha);
npoly = [1 -3 3+GG^3 -1];
n_roots = roots(npoly);
if isreal(n_roots(3)) && n_roots(3) < 1 && n_roots(3)>0
	n = n_roots(3);
elseif isreal(n_roots(2)) && n_roots(2) < 1 && n_roots(2)>0
	n = n_roots(2);
elseif isreal(n_roots(1)) && n_roots(1) < 1 && n_roots(2)>0
	n = n_roots(1);
else
	n_roots;
	n = 0.3;
	warning('cubic eqn has no good solution, force nbar.');
end
y = A*k^aalpha*n^(1-aalpha);
v = ttheta*(e+y);
m = v*(e+y)/(v^iiota+(e+y)^iiota)^(1/iiota);

% Find mk and me;
mk = (1-ddelta+U*aalpha*y/k)/c;
me = U/c;

% Find next state
eplus = e + y - m;
I = y-c+e-eplus-kkappa_S*v-kkappa_F*m-eplus*h-gyratio*y;
kplus = (1-ddelta)*k + I;

% pack stuff into control
control.kplus = kplus;
control.eplus = eplus;
control.mk = mk;
control.me = me;
control.q = q;
control.n = n;

end