function residual = solve_ttheta(ttheta,param)
ttau = param.ttau;
kkappa_F = param.kkappa_F;
kkappa_S = param.kkappa_S;
iiota = param.iiota;
bbeta = param.bbeta;
h = param.h;

f = @(ttheta) (1+ttheta^(-iiota))^(-1/iiota);
q = @(ttheta) (1+ttheta^(iiota))^(-1/iiota);

residual = (1-kkappa_F/ttau-kkappa_S/ttau/q(ttheta))/bbeta ...
	-( f(ttheta)*(1-kkappa_F)-kkappa_S*ttheta - h )/(1-bbeta*(1-f(ttheta)));
end