clear; clc; format long

%% preparation

% d      = 2;          % number of dimensions
% pi_pdf = repmat(ERADist('standardnormal','PAR'), d, 1);   % n independent rv

%% limit state function

a = 0.05;  b = 0.18;

g  = @(x) 5*(4-2.1*(a*x(:,1)).^2 + (a*x(:,1)).^4./3).*(a*x(:,1)).^2 + 5*(a*x(:,1)).*(b*x(:,2)) +10*((b*x(:,2)).^2-1).*(b*x(:,2)).^2 + 2.6;

%% MCS

% mu = zeros(1,d); stdu = eye(d);         % mean and standard deviation of input variable
% n_mcs = 10^6;
% 
% for i = 1:1000
%    u = lhsnorm(mu,stdu,n_mcs);
%    y = g(u);
%    pf_mcs(i) = sum(y <= 0)/n_mcs;
% end
% 
% pf_ref = mean(pf_mcs);  
 
pf_ref = 3.713900000000000e-05

%% Sequential directional importance sampling

nf     = 150;  % importance directions per level 
len    = 5;    % length of each Markov chain 
sigma  = 3;    % initial sigma
tarCoV = 1.5;  % target coefficient of variation of important weight
num    = 1;  % number of runs

for i = 1 : num                                                           % repeated runs
    i
   [pf(i), pf1(i), cov(i), n_cost(i),level_SuS(i),level_SDIS(i),cov_SuS(i),cov_SDIS(i)] = SDIS_figure(g,pi_pdf,nf,len,sigma,d,tarCoV);  % run SDIS algorithm
end

n_m  = mean(n_cost')        % mean of computational costs
pf_m = mean(pf')            % mean of failure probability
cv_m = mean(cov')           % mean of coefficient of variation
cv   = std(pf')./mean(pf')  % coefficient of variation of multiple runs

mse = (pf_m-pf_ref)^2 + var(pf);
eff = pf_ref*(1-pf_ref)/mse/n_m;


% corr = [0.8 0.4 0.2 0.1 0.05];
% rou   = 0;
% for i = 1:5
%     rou = rou + 2*(1-i/10^4)*corr(i);
% end

