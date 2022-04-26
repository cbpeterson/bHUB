
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This DEMO illustrates Bayeisan inference of hub nodes across multiple networks
%%% An illustrative example in Kim et al (2018+)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
addpath './Wang';
	
%% Data for an illustrative example from Kim et al (2018+)	
% K is number of subgroups
K = 3;	

% X1,X2 and X3 are data of the 60-node
X1 = textread('./X1.txt');
X2 = textread('./X2.txt');
X3 = textread('./X3.txt');

% p is number of variables
p = size(X1, 2);

% n is number of samples
n = size(X1, 1);	

% X'X matrix (p x p matrix)
S1 = X1' * X1;
S2 = X2' * X2;
S3 = X3' * X3;


% normalize matrix 
S1a = S1;
S2a = S2;
S3a = S3;
for i = 1:p
	for j = 1:p
		S1a(i,j) = S1(i,j)/sqrt(S1(i,i)*S1(j,j));
		S2a(i,j) = S2(i,j)/sqrt(S2(i,i)*S2(j,j));
		S3a(i,j) = S3(i,j)/sqrt(S3(i,i)*S3(j,j));
	end
end
S1 = S1a * n;
S2 = S2a * n;	
S3 = S3a * n;

	
% Number of MCMC iterations before burnin
burnin = 2000; 
	
% Number of MCMC iterations after burnin
nmc = 10000;

	
% Initial values for covariance matrix for each group
SS = eye(p);


%% Apply bHUB
[H_save,C_save,Z_save,ar_hub,ar_tau,ar_b2] = ...
			hubGGMMCMC(cat(3, S1, S2, S3), repmat(n, 1, K), repmat(SS, [1, 1, K]), burnin,nmc);

% Compute the posterior probability of inclusion for each edge, and the posterior
% probability of each node being a hub
ppi_edges = mean(Z_save, 4);
ppi_hub = mean(H_save, 3);


%% Print results
dlmwrite('./ppi_hub.txt', ppi_hub, '\t')
for g = 1:K
	fname =  sprintf('./ppi_edgesGroup%d.txt',g);
	dlmwrite(fname, ppi_edges(:,:,g), '\t')
end
