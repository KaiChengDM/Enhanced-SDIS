clear; clc; format long

%% preparation

d      = 2;          % number of dimensions
pi_pdf = repmat(ERADist('standardnormal','PAR'), d, 1);   % n independent rv

%% limit state function

g  = @(x)30./((4*(x(:,1)+2).^2/9+x(:,2).^2/25).^2 + 1) + 20./(((x(:,1)-2.5).^2/4 + (x(:,2)-0.5).^2/25).^2+1) - 5

%% MCS

% mu = zeros(1,d); stdu = eye(d);         % mean and standard deviation of input variable
% n_mcs = 10^6;
% 
% for i = 1:100
%    u = lhsnorm(mu,stdu,n_mcs);
%    y = g(u);
%    pf_mcs(i) = sum(y <= 0)/n_mcs;
% end
% 
% pf_ref = mean(pf_mcs);  

pf_ref = 1.125000000000000e-05

%% Sequential directional importance sampling

nf     = 150;  % importance directions per level 
len    = 5;    % length of each Markov chain 
sigma  = 3;    % initial sigma
tarCoV = 1.5;  % target coefficient of variation of important weight
num    = 300;  % number of runs

for i = 270: num                                                           % repeated runs
    i
   [pf(i), pf1(i), cov(i), n_cost(i),level_SuS(i),level_SDIS(i),cov_SuS(i),cov_SDIS(i)] = SDIS(g,pi_pdf,nf,len,sigma,d,tarCoV);  % run SDIS algorithm
end

n_m  = mean(n_cost')        % mean of computational costs
pf_m = mean(pf')            % mean of failure probability
cv_m = mean(cov')           % mean of coefficient of variation
cv   = std(pf')./mean(pf')  % coefficient of variation of multiple runs

mse = (pf_m-pf_ref)^2 + var(pf);
eff = pf_ref*(1-pf_ref)/mse/n_m;

