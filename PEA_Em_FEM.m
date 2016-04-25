%% Housekeeping
clear
close all
clc
format long
addpath(genpath('./tools'))
addpath(genpath('./param'))

%% Set the stage
mypara;
nA = 9;
nK = 50;
nE = 50;
T = 10000;
[P,lnAgrid] = rouwen(rrho_z,0,ssigma_z/sqrt(1-rrho_z^2),nA);
Anodes = exp(lnAgrid);
P = P';
min_lnA = lnAgrid(1); max_lnA = lnAgrid(end);
min_K = 0.01; max_K = 50;
min_E = 0.01; max_E = 40;
damp_factor = 0.1;
maxiter = 10000;
tol = 1e-6;
options = optimoptions(@fsolve,'Display','none','Jacobian','off');

%% Grid creaton
ln_kgrid = linspace(log(min_K),log(max_K),nK);
Knodes = exp(ln_kgrid);
ln_egrid = linspace(log(min_E),log(max_E),nE);
Enodes = exp(ln_egrid);
N = nA*nK*nE;
[Kmesh,Amesh,Nmesh] = meshgrid(Knodes,Anodes,Enodes);

% %% Encapsulate all parameters
% param = [... 
%  bbeta; % 1
%  ggamma; % 2
%  kkappa; % 3
%  eeta; % 4
%  rrho; %5
%  ssigma; %6
%  x; % 7
%  aalpha; % 8
%  ddelta; % 9
%  xxi; % 10
%  ttau; % 11
%  z % 12
%  ];

%% Precomputation and initial guess
tot_stuff = zeros(N,1); ustuff = zeros(N,1);
EMKval = zeros(nA,nK,nE); EMEval = EMKval;
EMKval_temp = EMKval; EMEval_temp = EMEval;
% parfor i = 1:N
% 	[i_a,i_k,i_n] = ind2sub([nA,nK,nE],i);
% 	a = Anodes(i_a); k  = Knodes(i_k); n = Nnodes(i_n); %#ok<PFBNS>
% 	tot_stuff(i) = a*k^aalpha*n^(1-aalpha) + (1-ddelta)*k + z*(1-n);
% 	ustuff(i) = xxi*(1-n)^(1-eeta);
% end
if (exist('PEA_Em_FEM.mat','file'))
	load('PEA_Em_FEM.mat','EMKval','EMEval');
else
	EMKval = zeros(nA,nK,nE); EMEval = EMKval;
	EMKval_temp = EMKval; EMEval_temp = EMEval;
	coeff_lnmk = zeros(4,1); coeff_lnme = zeros(4,1);
	coeff_lnmk(1) = 1.9533;
	coeff_lnmk(2) = -0.2407;
	coeff_lnmk(3) = -0.5138;
	coeff_lnmk(4) = -0.0574;
	
	coeff_lnme(1) = 1.8187;
	coeff_lnme(2) = -0.2634;
	coeff_lnme(3) = -0.4928;
	coeff_lnme(4) = -0.0668;
	parfor i = 1:N
		[i_a,i_k,i_e] = ind2sub([nA,nK,nE],i);
		a = Anodes(i_a); k  = Knodes(i_k); e = Enodes(i_e); %#ok<PFBNS>
		EMKval(i) = exp([1 log(a) log(k) log(e)]*coeff_lnmk);
		EMEval(i) = exp([1 log(a) log(k) log(e)]*coeff_lnme);
	end
end
	



%% Iteration
diff = 10; iter = 0;
while (diff>tol && iter <= maxiter)
	% Pack grids, very important
	grids.EMKval = EMKval;
	grids.EMEval = EMEval;
	grids.Knodes = Knodes;
	grids.Enodes = Enodes;
	
    %% Time iter step, uses endo grid technique
    parfor i = 1:N
		
        [i_a,i_k,i_e] = ind2sub([nA,nK,nE],i);
        e = Enodes(i_e); k = Knodes(i_k); A = Anodes(i_a);
		state = [A k e];
		
		% Find current control vars
		control = state2control_FEM(state,i_a,grids,param);
		kplus = control.kplus;
		eplus = control.eplus;
		
		% Find the expected EM
        EMK_hat = 0; EME_hat = 0;
        for i_node = 1:nA
            aplus = Anodes(i_node);
			stateplus = [aplus kplus eplus];
			control_plus = state2control_FEM(stateplus,i_node,grids,param);
			
			EMK_hat = EMK_hat + P(i_a,i_node)*control_plus.mk;
			EME_hat = EME_hat + P(i_a,i_node)*control_plus.me;
        end
        
        EMKval_temp(i) = EMK_hat;
        EMEval_temp(i) = EME_hat;
    end
    
    %% Damped update
    EMKval_new = (damp_factor)*EMKval_temp+(1-damp_factor)*EMKval;
    EMEval_new = (damp_factor)*EMEval_temp+(1-damp_factor)*EMEval;
    
    %% Compute norm
    diff = norm([EMKval(:);EMEval(:)]-[EMKval_new(:);EMEval_new(:)],Inf);
    
    %% Update
    EMKval = EMKval_new;
    EMEval = EMEval_new;
    iter = iter+1;
    %% Display something
    iter
    diff

end;

%% Euler equation error
nk_ee = 60; nnn_ee = 60;
Kgrid = linspace(0.5*k_ss,1.5*k_ss,nk_ee);
Agrid = exp(lnAgrid);
Ngrid = linspace(0.96*n_ss,1.04*n_ss,nnn_ee);
EEerror_c = 999999*ones(nA,nk_ee,nnn_ee);
EEerror_v = 999999*ones(nA,nk_ee,nnn_ee);
cc = zeros(nA,nk_ee,nnn_ee);
vv = zeros(nA,nk_ee,nnn_ee);
tthetattheta = zeros(nA,nk_ee,nnn_ee);
cc_dynare = cc;
vv_dynare = vv;
tthetattheta_dynare = tthetattheta; 

for i_a = 1:nA
    a = Agrid(i_a);
    for i_k = 1:nk_ee
        k = Kgrid(i_k);
        for i_e = 1:nnn_ee
            e = Ngrid(i_e);
			tot_stuff = a*k^aalpha*e^(1-aalpha)+(1-ddelta)*k+z*(1-e);
			ustuff = xxi*(1-e)^(1-eeta);
            
            EMK = globaleval(k,e,Knodes,Enodes,squeeze(EMKval(i_a,:,:)));
            EMF = globaleval(k,e,Knodes,Enodes,squeeze(EMEval(i_a,:,:)));
            c = 1/(bbeta*EMK);
            q = kkappa/c/(bbeta*EMF);            
            
            if q <= 0
                warning('q <= 0!!')
                q = 0;
                ttheta = 0;
                v = 0;
                kplus = tot_stuff - c - kkappa*v;
                nplus = (1-x)*e;
            else
                ttheta = (q/xxi)^(1/(eeta-1));
                v = ttheta*(1-e);
                kplus = tot_stuff - c - kkappa*v;
                nplus = (1-x)*e + xxi*v^eeta*(1-e)^(1-eeta);
            end
            
            cc(i_a,i_k,i_e) = c;
            cc_dynare(i_a,i_k,i_e) = exp( 2.130385+0.039519*(log(a)/rrho-0)+0.606879*(log(k)-log(k_ss))+0.005573*(log(e)-log(n_ss)) );
            vv(i_a,i_k,i_e) = v;
            vv_dynare(i_a,i_k,i_e) = exp( -2.899249+3.417972*(log(a)/rrho-0)+0.451375*(log(k)-log(k_ss))+(-17.928147)*(log(e)-log(n_ss)) );
            tthetattheta(i_a,i_k,i_e) = ttheta;
            tthetattheta_dynare(i_a,i_k,i_e) = exp( 0+3.417972*(log(a)/rrho-0)+0.451375*(log(k)-log(k_ss))+(-0.767653)*(log(e)-log(n_ss)) );

			% Find expected mh, mf tomorrow if current coeff applies tomorrow
            EMH_hat = 0; EME_hat = 0;
            for i_node = 1:nA
                aplus = Anodes(i_node);
                EMH_plus = globaleval(kplus,nplus,Knodes,Enodes,squeeze(EMKval(i_node,:,:)));
                EMF_plus = globaleval(kplus,nplus,Knodes,Enodes,squeeze(EMEval(i_node,:,:)));
                cplus = 1/(bbeta*EMH_plus);
                qplus = kkappa/cplus/(bbeta*EMF_plus);
				if qplus <= 0
					% warning('qplus <= 0!!')
					qplus = 0;
					tthetaplus = 0;
					vplus = 0;
					EMH_hat = EMH_hat + P(i_a,i_node)*((1-ddelta+aalpha*aplus*(kplus/nplus)^(aalpha-1))/cplus);
					EME_hat = EME_hat + P(i_a,i_node)*(( (1-ttau)*((1-aalpha)*aplus*(kplus/nplus)^aalpha-z-ggamma*cplus) )/cplus );
				else
					tthetaplus = (qplus/xxi)^(1/(eeta-1));
					vplus = tthetaplus*(1-nplus);
					EMH_hat = EMH_hat + P(i_a,i_node)*((1-ddelta+aalpha*aplus*(kplus/nplus)^(aalpha-1))/cplus);
					EME_hat = EME_hat + P(i_a,i_node)*(( (1-ttau)*((1-aalpha)*aplus*(kplus/nplus)^aalpha-z-ggamma*cplus) + (1-x)*kkappa/qplus - ttau*kkappa*tthetaplus )/cplus );
				end
            end

			c_imp = 1/(bbeta*EMH_hat);
			q_imp = kkappa/c_imp/(bbeta*EME_hat);
			ttheta_imp = (q_imp/xxi)^(1/(eeta-1));
			v_imp = ttheta_imp*(1-e);

            EEerror_c(i_a,i_k,i_e) = abs((c-c_imp)/c_imp);   
            EEerror_v(i_a,i_k,i_e) = abs((v-v_imp)/v_imp);  
        end
    end
end
EEerror_c_inf = norm(EEerror_c(:),inf)
EEerror_v_inf = norm(EEerror_v(:),inf)

EEerror_c_mean = mean(EEerror_c(:));
EEerror_v_mean = mean(EEerror_v(:));

%% Simulation
P_cdf = cumsum(P,2);
aindexsim = zeros(1,T); aindexsim(1) = ceil(nA/2);
ksim = kbar*ones(1,T); esim = ebar*ones(1,T);
% tthetasim = zeros(1,T); vsim = zeros(1,T); usim = zeros(1,T);
for t = 1:T
    asim(t) = Anodes(aindexsim(t)); a = asim(t);
    k = ksim(t); e = nsim(t);
	state = [a k e];
	
	control = state2control_FEM(state,aindexsim(t),grids,param);
	
    tot_stuff = a*k^aalpha*e^(1-aalpha)+(1-ddelta)*k+z*(1-e);
    ustuff = xxi*(1-e)^(1-eeta);
    
    EMK = globaleval(k,e,Knodes,Enodes,squeeze(EMKval(aindexsim(t),:,:)));
    EMF = globaleval(k,e,Knodes,Enodes,squeeze(EMEval(aindexsim(t),:,:)));
    c = 1/(bbeta*EMK);
    q = kkappa/c/(bbeta*EMF);
    
    if q <= 0
        warning('q <= 0!!')
        q = 0;
        ttheta = 0;
        v = 0;
        kplus = tot_stuff - c - kkappa*v;
        nplus = (1-x)*e;
    else
        ttheta = (q/xxi)^(1/(eeta-1));
        v = ttheta*(1-e);
        kplus = tot_stuff - c - kkappa*v;
        nplus = (1-x)*e + xxi*v^eeta*(1-e)^(1-eeta);
    end
    
    tthetasim(t) = ttheta;
    vsim(t) = v;
    usim(t) = 1-nsim(t);
    
    if t <= T-1
        uu = rand;
        aindexsim(t+1) = find(P_cdf(aindexsim(t),:)>=uu,1,'first');
        ksim(t+1) = kplus;
        nsim(t+1) = nplus;
    end
end

[~,ttheta_cyc] = hpfilter(log(tthetasim),12^2*100);
[~,u_cyc] = hpfilter(log(usim),12^2*100);
[~,v_cyc] = hpfilter(log(vsim),12^2*100);

%% Export results
mkdir('results')
h_c = figure;
plot(Kgrid,squeeze(EEerror_c(ceil(nA/2),:,ceil(nnn_ee/2))))
title('Euler Error of Consumption')
print(h_c,'-dpsc','./results/EEerror_c.eps')

h_v = figure;
plot(Kgrid,squeeze(EEerror_v(ceil(nA/2),:,ceil(nnn_ee/2))))
title('Euler Error of Vacancy')
print(h_v,'-dpsc','./results/EEerror_v.eps')

result_mf = @(k,n) globaleval(k,n,Knodes,Enodes,squeeze(EMEval(1,:,:)));
h_EMF = figure;
ezsurf(result_mf,[Kgrid(1),Kgrid(end),Ngrid(1),Ngrid(end)])
print(h_EMF,'-dpsc','./results/EMF.eps')

v_policy = figure;
plot(Ngrid,squeeze(vv(1,1,:)))
title('Vacancy policy at lowerest productivity and capital.')
print(v_policy,'-dpsc','./results/v_policy.eps')

c_policy = figure;
plot(Kgrid,squeeze(cc(ceil(nA/2),:,ceil(nnn_ee/2))),Kgrid,squeeze(cc_dynare(ceil(nA/2),:,ceil(nnn_ee/2))))
title('Consumption policies at SS.')
print(c_policy,'-dpsc','./results/c_policy.eps')
xlabel('Capital')

ttheta_policy = figure;
plot(Kgrid,squeeze(tthetattheta(ceil(nA/2),:,ceil(nnn_ee/2))),Kgrid,squeeze(tthetattheta_dynare(ceil(nA/2),:,ceil(nnn_ee/2))))
title('\theta around SS')
print(ttheta_policy,'-dpsc','./results/ttheta_policy.eps')
xlabel('Capital')

ttheta_policyN = figure;
plot(Ngrid,squeeze(tthetattheta(ceil(nA/2),ceil(nk_ee/2),:)),Ngrid,squeeze(tthetattheta_dynare(ceil(nA/2),ceil(nk_ee/2),:)))
title('\theta around SS')
print(ttheta_policyN,'-dpsc','./results/ttheta_policy2.eps')
xlabel('Employment')

ttheta_policyA = figure;
plot(Anodes,squeeze(tthetattheta(:,ceil(nk_ee/2),ceil(nnn_ee/2))),Anodes,squeeze(tthetattheta_dynare(:,ceil(nk_ee/2),ceil(nnn_ee/2))))
title('\theta around SS')
xlabel('Productivity')
print(ttheta_policyA,'-dpsc','./results/ttheta_policy3.eps')

save('PEA_Em_FEM.mat');
