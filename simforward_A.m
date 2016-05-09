function path = simforward_A(init_state,Ainit_idx,innov,periods,grids,param)
%%load stuff
path = zeros(5,1+periods); % now we only focus on CIPI
rrho_z = param.rrho_z;
ssigma_z = param.ssigma_z;
Anodes = grids.Anodes;
P_cdf = grids.P_cdf;
aalpha = param.aalpha;

%% use minus one period state to find period one states
control_minusone = state2control_FEM_simple(init_state,Ainit_idx,grids,param);
ln_Aone = rrho_z*log(init_state(1))+ssigma_z*innov;
Aone = exp(ln_Aone);
[~,i_Aone] = min(abs(Aone-Anodes));
state(1) = Aone;
state(2) = control_minusone.kplus;
state(3) = control_minusone.eplus;
A_idx = i_Aone;

%% Compute initial IRF variable as benchmark
Y_benchmark = init_state(1)*init_state(2)^aalpha;
CIPI_benchmark = state(3) - init_state(3);

%% simulate forward
for t = 1:length(path)
	% Find control variales and IRF of interest
	control = state2control_FEM_simple(state,A_idx,grids,param);
	CIPI = control.eplus - state(3);
	q = control.q;
	f = control.f;
	v = control.v;
	Y = state(1)*state(2)^aalpha;
	path(1,t) = CIPI;
	path(2,t) = Y;
	path(3,t) = CIPI/Y;
	path(4,t) = q;
	path(5,t) = f;
	path(6,t) = v;
	
	% update future state
	uu = rand;
    A_idx = find(P_cdf(A_idx,:)>=uu,1,'first');
	state(1) = Anodes(A_idx);
	state(2) = control.kplus; state(3) = control.eplus;
	
end

end