function [H_save,C_save,Z_save,ar_hub,ar_tau,ar_b2] = hubGGMMCMC(S,n,Sig,burnin,nmc)
% Input parameters:
%   S: p x p x K array of cross product matrix:  Y*Y' if Y is a p x n matrix of random sample
%   n: vector of sample sizes for each group
%   Sig: p x p x K array of initial guess of Sigma
%   burnin: number of burn-in iterations
%   nmc: number of MCMC sample iterations saved
% Output parameters:
%   H_save: p x K x nmc sample of hub-specific parameters
%   C_save: p x p x K x nmc sample of precision matrix
%   Z_save: p x p x K x nmc sample of graph matrix
%   ar_hub: p x K matrix of cceptance rate for hub indicators
%   ar_tau: a scalar of acceptance rate for tau
%   ar_b2: a scalar of acceptance rate for beta2


	% K is number of sample groups
	K = size(n, 2);

	% p is number of variables
	p = size(S, 1);

	% Set up matrices for return values
	% hub indicator
	H_save = zeros(p, K, nmc);

	% precision matrices
	C_save = zeros(p, p, K, nmc);

	% graph matrices
	Z_save = zeros(p, p, K, nmc);

	% tau 
	tau_save = zeros(2, 1, nmc);
 
	% beta2 
	b2_save = zeros(1, nmc);
 
	% Acceptance rate for gamma 
	ar_tau = zeros(1, 1);

	% Acceptance rate for each hub indicator
	ar_hub = zeros(p, K);

	% Acceptance rate for beta2
	ar_b2 = 0;

	% C: p x p x K array of initial precision matrices for each group
	C = zeros(p, p, K);
	% Initial graphs matrices for each group
	Z = zeros(p, p, K);
	% Initial hub indicators for each group
	H = zeros(p, K);


	% Prior for variances
	% 	V0:  p x p matrix of the  small variance  components;  V0(i,j) is the small variance of the precision element C(i,j)
	% 	V1:  p x p matrix of the  large variance  components;  V1(i,j) is the small variance of the precision element C(i,j)
	% 	lambda: hyperparameter for the diagonal element Usually set to be 1;

	v0 = 0.03;
	h = 100;
	v1 = h*v0;
    lambda = 1;
    V0 = v0*ones(p);
    V1 = v1*ones(p);

	% eta: baseline sparsity for each node's hubness
	% logit^(-1)(eta) = 0.01
	eta = log(1/99); 
	
	% Prior for edge inclusion
	%   beta1: prior marginal edge "inclusion probability" for small probability
	%   beta2: prior marginal edge "inclusion probability" for large probability
	beta1 = 0.001;
	beta2 = 0.001;

	% mean(beta2) = 0.025
	a_prop = 152.319;
	b_prop =  5940.441;

	% Initial prior for tau 
	tau = 0;
	tau_prop = 0;


	% Proposal parameters for MH steps (tau)
	a_tau = 1.5; 
	b_tau = 0.2; 
	tau_sig = 0.0005; %variance of truncated normal proposal

	% Perform MCMC sampling
	for iter = 1: burnin+nmc    
            
		if(mod(iter, 1000) == 0)
			fprintf('iter = %d\n', iter);
		end
	
    	% Update graph and precision matrix for each group using modified code from Hao Wang
    	for cur_graph = 1:K

			% Updating precision matrix and graph for current graph
			[C(:,:,cur_graph),Z(:,:,cur_graph),Sig(:,:,cur_graph)] = hubGGM_SSVS(S(:,:,cur_graph),n(cur_graph),Sig(:,:,cur_graph),V0,V1,1,H(:,cur_graph),beta1,beta2,0,1);
		end	
	
		% Updating node indicators 
		for k = 1:K
			no_k = setdiff(1:K, k);

			for i = 1:p
				ind_noi = setdiff(1:p,i);

				% Calculate MH ratio
				H1 = sum(H(i,no_k));
				h_j = H(ind_noi, k);
				
				A = eta + 2 * tau * H1; 

				pii = beta1*(1-h_j) + beta2*h_j;

				hub_one = beta2^sum(Z(i,ind_noi,k)) * (1-beta2)^sum(1 - Z(i,ind_noi,k)) * exp(A)/(1+exp(A)); 
				hub_zero = prod(pii'.^ Z(i,ind_noi,k)) * prod((1-pii)'.^ (1-Z(i,ind_noi,k))) * 1/(1+exp(A));
	        	H(i, k) = binornd(1, hub_one/(hub_one + hub_zero));
			end
		end

		% Updating beta2
		beta_a = 0;
		beta_b = 0;
		for k = 1:K 
			for i = 1:(p-1)
				for j = (i+1):p
					temp = (1 - H(i, k)) * (1 - H(j, k)); 
					beta_a = beta_a + (1 - temp) * Z(i,j,k);
					beta_b = beta_b + (1 - temp) * (1-Z(i,j,k));
				end
			end
		end

        % sample from beta posterior
		beta2 = betarnd(a_prop + beta_a, b_prop + beta_b);

	
		% Updating tau's
		tau_prop = gamrnd(a_tau,b_tau);
		M = ones(K,K) - eye(K); 
		lpost = 0;
		for i = 1:p
			lpost = lpost + ...
			log(calc_mrf_C(M, eta, tau)) - log(calc_mrf_C(M, eta, tau_prop)) + ...
			(tau_prop - tau) .* H(i, :) * M * transpose(H(i, :));
		end
		
		% Ga(a_tau,b_tau) and truncated normal walk
		logd = (a_tau - 1).*(log(tau_prop)-log(tau)) - (tau_prop-tau)./b_tau;
		logwalk = (tau_prop - tau)^2./(2 * tau_sig);
		% calculating MH ratio
		if (log(unifrnd(0,1,1,1)) < (lpost + logd + logwalk))
			tau = tau_prop;
		 	ar_tau = ar_tau + 1 / (burnin + nmc);
		end

        % Save current values for return if we are past the burn-in
		if iter > burnin
			H_save(:,:,iter-burnin) = H(:,:);
			C_save(:,:,:,iter-burnin) = C(:,:,:);
			Z_save(:,:,:,iter-burnin) = Z(:,:,:);
		end
	end

	
end	
