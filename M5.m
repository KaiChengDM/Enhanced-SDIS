
clear;  clc;

%% limit state function

dim = [10 100 1000];          % number of dimensions

for k = 1:3
    k
    d = dim(k);                   % number of dimensions
   
    pi_pdf = repmat(ERADist('standardnormal','PAR'), d, 1);   % n independent rv

    Beta = 4;                     % reliability index
    g1   = @(x)Beta*sqrt(d) - sum(x');
    g2   = @(x)Beta*sqrt(d) + sum(x');
    g    = @(x)min(g1(x),g2(x));  % limit state function  
    pf_ref = 2*normcdf(-Beta)     % true failure probability

%% Sequential directional importance sampling

nf     = 150;  % importance directions per level 
len    = 5;    % length of each Markov chain 
sigma  = 3;    % initial sigma
tarCoV = 1.5;  % target coefficient of variation of important weight
num    = 300;  % number of runs

for i = 1: num                                                           % repeated runs
    i
   [pf(i), pf1(i), cov(i), n_cost(i),level_SuS(i),level_SDIS(i),cov_SuS(i),cov_SDIS(i)] = SDIS(g,pi_pdf,nf,len,sigma,d,tarCoV);  % run SDIS algorithm
end

   n_m(k)   = mean(n_cost')        % mean of computational costs
   pf_m(k)  = mean(pf')            % mean of failure probability
   pf_m1(k) = mean(pf1')          % mean of failure probability
   cv_m(k)  = mean(cov')           % mean of coefficient of variation
   cv(k)    = std(pf')./mean(pf')  % coefficient of variation of multiple runs

   L_SuS(k)  = mean(level_SuS);
   L_SDIS(k) = mean(level_SDIS);

   mse(k) = (pf_m(k)-pf_ref)^2 + var(pf);
   eff(k) = pf_ref*(1-pf_ref)/mse(k)/n_m(k);

end
