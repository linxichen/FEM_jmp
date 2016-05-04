function control = state2control_simple(state,coeff_mk,coeff_me,param)
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

% find expectation terms
EM = exp([1 log(state)]*[coeff_mk coeff_me]);
A = state(1); k = state(2); e = state(3);

% Find expectation terms tomorrow
Eme = EM(2); Emk = EM(1);

% Back out other variables
c = (bbeta*Emk)^(-1);
q = kkappa_S/ttau*(1-c*bbeta*Eme-kkappa_F/ttau)^(-1);
Pm = 1 - kkappa_F - kkappa_S/q;
if q >= 1
	q = 1.0-1e-7;
	% warning('q>=1');
elseif q <=0
	q = 1e-7;
	% warning('q<=0');
end
if q ~= 0
	ttheta = (q^(-iiota) - 1)^(1/iiota);
else
	ttheta = 1e7;
end
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
eplus = max(0.1,eplus);
I = y-c+e-eplus-kkappa_S*v-kkappa_F*m-e*h;
kplus = (1-ddelta)*k + I;
kplus = max(0.1,kplus);

% pack stuff into control
control.kplus = kplus;
control.eplus = eplus;
control.mk = mk;
control.me = me;
control.q = q;
control.f = f;
control.n = 1.0;
control.CIPI = eplus - e;

end