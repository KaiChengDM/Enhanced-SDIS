
clear;  clc; format long

%% Preparation

dim = [10 100 1000];          % number of dimensions

pf_const = 5*10^-5;

niu  = 1;

%% MCS

% mu = zeros(1,d); stdu = eye(d);         % mean and standard deviation of input variable
% n_mcs = 10^6;
% 
% for i = 1:1
%    u = lhsnorm(mu,stdu,n_mcs);
%    y = g(3*u);
%    pf_mcs(i) = sum(y <= 0)/n_mcs;
% end

%% Sequential directional importance sampling

 nf     = 150;      % Importance directions per level 
 len    = 5;        % Length of each Markov chain 
 sigma  = 3;        % Initial sigma
 tarCoV = 1.5;      % target coefficient of variation of important weight

for k = 2:3

    d   = dim(k);

    pi_pdf = repmat(ERADist('standardnormal','PAR'), d, 1);   % n independent rv

    C_a = gaminv(1-pf_const,d,niu);
    g0  = @(x)C_a + sum(x');                 % Performance function
    g   = @(x)g0(niu.*log(normcdf(-x)));

    pf_ref = 1 - gamcdf(C_a,d,niu)           % True failure probability

    for i = 1:300
        i
        [pf(i), pf1(i),cov(i), n_cost(i),level_SuS(i),level_SDIS(i),cov_SuS(i),cov_SDIS(i)] = SDIS(g,pi_pdf,nf,len,sigma,d,tarCoV);
 
    end

    N(k)   = mean(n_cost');
    Pf(k)  = mean(pf');
    Pf1(k) = mean(pf1');
    CV1(k) = mean(cov');
    CV2(k) = std(pf')./mean(pf');

    L_SuS(k)  = mean(level_SuS);
    L_SDIS(k) = mean(level_SDIS);

    mse(k) = (Pf(k)-pf_ref)^2 + var(pf);
    eff(k) = pf_ref*(1-pf_ref)/mse(k)/N(k);

end


