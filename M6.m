clear; clc;  

%% preparation

d      = 2;          % number of dimensions
pi_pdf = repmat(ERADist('standardnormal','PAR'), d, 1);   % n independent rv

%% limit state function

Wp = @(x)sqrt(x(:,3)./x(:,1));
Ws = @(x)sqrt(x(:,4)./x(:,2));
Wa = @(x)(Wp(x)+Ws(x))./2;
Fa = @(x)(x(:,5)+x(:,6))./2;
gm = @(x)x(:,2)./x(:,1); 
St = @(x)(Wp(x)-Ws(x))./Wa(x);
g0 = @(x)x(:,7)-3.*x(:,4).*sqrt(pi.*x(:,8)./(4.*x(:,6).*Ws(x).^3).*(Fa(x).*x(:,6)./(x(:,5).*x(:,6).*(4*Fa(x).^2+St(x).^2)+gm(x).*Fa(x).^2)).*(x(:,5).*Wp(x).^3+x(:,6).*Ws(x).^3).*Wp(x)./(4.*Fa(x).*Wa(x).^4));

mu_Input  = [1.5 0.01 1 0.01 0.05 0.02 21.5 100];
cov_Input = [0.1 0.1 0.2 0.2 0.4 0.5 0.1 0.1];

sigma_Input = mu_Input.*cov_Input;  d = 8;

mean_ln  = log(mu_Input)-0.5.*log(1+sigma_Input.^2./mu_Input.^2); % 4.42e-5
sigma_ln = sqrt(log(1+sigma_Input.^2./mu_Input.^2));    % SUS:4.3803e-05, 0.4663,5.1620e+03

g = @(x)g0(exp(x.*sigma_ln + mean_ln));  % Transform performance function into standard normal  space

pf_ref = 4.42e-5;

%% MCS

% mu    = zeros(d,1);  
% sigma = eye(d); 
% n_mcs = 10^6;
% 
% for i = 1:1000
%     i
%    u = lhsnorm(mu,sigma,n_mcs);
%    y = g(u);
%    pf_mcs(i) = sum(y <= 0)/n_mcs;
% end
% 
% pf_ref = mean(pf_mcs); 

%% Sequential directional importance sampling

nf     = 150;  % importance directions per level 
len    = 5;    % length of each Markov chain 
sigma  = 3;    % initial sigma
tarCoV = 1.5;  % target coefficient of variation of important weight
num    = 300;  % number of runs

for i = 1: num  % repeated runs
    i
   [pf(i), pf1(i), cov(i), n_cost(i),level_SuS(i),level_SDIS(i),cov_SuS(i),cov_SDIS(i)] = SDIS(g,pi_pdf,nf,len,sigma,d,tarCoV);  % run SDIS algorithm
end

n_m  = mean(n_cost')        % mean of computational costs
pf_m = mean(pf')            % failure probability
cv_m = mean(cov')           % mean of coefficient of variation
cv   = std(pf')./mean(pf')  % coefficient of variation of multiple runs

mse = (pf_m-pf_ref)^2 + var(pf);
eff = pf_ref*(1-pf_ref)/mse/n_m

