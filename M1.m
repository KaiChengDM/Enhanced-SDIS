
clear; close all; clc; format long
% rng(10)

%% definition of the random variables

d      = 2;          % number of dimensions
pi_pdf = repmat(ERADist('standardnormal','PAR'), d, 1);   % n independent rv

%% limit state function  

g1 = @(x)3+0.1.*(x(:,1)-x(:,2)).^2-(x(:,1)+x(:,2))./2^0.5;
g2 = @(x)3+0.1.*(x(:,1)-x(:,2)).^2+(x(:,1)+x(:,2))./2^0.5;
g3 = @(x)(x(:,1)-x(:,2))+6./2^0.5;
g4 = @(x)(x(:,2)-x(:,1))+6./2^0.5;
g  = @(x)min([g1(x),g2(x),g3(x),g4(x)]')'+2;  

%% Sequential directional importance sampling

nf     = 150;  % importance directions per level 
len    = 5;    % length of each Markov chain 
sigma  = 3;    % initial sigma
tarCoV = 1.5;  % target coefficient of variation of important weight
num    = 3;    % number of runs

for i = 1 : num                                                           % repeated runs
    i
   [pf(i), pf1(i), cov(i), n_cost(i),level_SuS(i),level_SDIS(i),cov_SuS(i),cov_SDIS(i)] = SDIS(g,pi_pdf,nf,len,sigma,d,tarCoV);  % run SDIS algorithm
end

n_m  = mean(n_cost')        % mean of computational costs
pf_m = mean(pf')            % mean of failure probability
cv_m = mean(cov')           % mean of coefficient of variation
cv   = std(pf')./mean(pf')  % coefficient of variation of multiple runs