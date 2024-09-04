clear; clc; format long

%% preparation

d      = 2;          % number of dimensions
pi_pdf = repmat(ERADist('standardnormal','PAR'), d, 1);   % n independent rv

%% limit state function
 
g_fun = @(x) min([3.2 + (1/sqrt(d))*(x(:,1)+x(:,2)), ...
                0.1*(x(:,1)-x(:,2)).^2 - (x(:,1)+x(:,2))./sqrt(d) + 2.5 ], [], 2)+3;  
g     = @(x)g_fun(x);    % limit state function

pf_ref = 1.10e-8;

%% Sequential directional importance sampling

nf     = 150;  % importance directions per level 
len    = 5;    % length of each Markov chain 
sigma  = 3;    % initial sigma
tarCoV = 1.5;  % target coefficient of variation of important weight
num    = 300;  % number of runs

for i = 1 : num                                                            % repeated runs
    i
   [pf(i), pf1(i), cov(i), n_cost(i),level_SuS(i),level_SDIS(i),cov_SuS(i),cov_SDIS(i)] = SDIS(g,pi_pdf,nf,len,sigma,d,tarCoV);  % run SDIS algorithm
end

n_m  = mean(level_SuS')        % mean of computational costs
pf_m = mean(pf')            % mean of failure probability
cv_m = mean(cov')           % mean of coefficient of variation
cv   = std(pf')./mean(pf')  % coefficient of variation of multiple runs

mse= (pf_m-pf_ref)^2 + var(pf);
eff = pf_ref*(1-pf_ref)/mse/n_m;  % compute the relative efficiency

