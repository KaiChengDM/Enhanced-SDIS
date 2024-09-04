
clear;  clc; format long

%% Limit state function

dim = [10 100 1000];          % number of dimensions

kappa = 1;
Beta  = 3.5;

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
% pf_ref = mean(pf_mcs);  % 1.369642857142857e-05

 pf_ref = 2.9173e-04;
 
%% Sequential directional importance sampling
 
 nf     = 150;      % Importance directions per level 
 len    = 5;        % Length of each Markov chain 
 sigma  = 3;        % Initial sigma
 tarCoV = 1.5;      % target coefficient of variation of important weight

for k = 1:3
   
    d = dim(k);
    pi_pdf = repmat(ERADist('standardnormal','PAR'), d, 1);   % n independent rv

    g1   = @(x)-sum(x')'/sqrt(d) + Beta + kappa*(x(:,1)-x(:,2)).^2/10;
    g2   = @(x) sum(x')'/sqrt(d) + Beta + kappa*(x(:,1)-x(:,2)).^2/10;
    g    = @(x)min(g1(x),g2(x));  

    for i = 1: 300
        i
        [pf(i), pf1(i),cov(i), n_cost(i), level_SuS(i),level_SDIS(i),cov_SuS(i),cov_SDIS(i)]= SDIS(g,pi_pdf,nf,len,sigma,d,tarCoV);
       
    end

    N(k)   = mean(n_cost');
    Pf(k)  = mean(pf');
    pf_1(k) = mean(pf1');           % mean of failure probability
    CV1(k) = mean(cov');
    CV2(k) = std(pf')./mean(pf');
    
    L_SuS(k)  = mean(level_SuS);
    L_SDIS(k) = mean(level_SDIS);

    mse(k) = (Pf(k)-pf_ref)^2 + var(pf);
    eff(k) = pf_ref*(1-pf_ref)/mse(k)/N(k);

end
  
