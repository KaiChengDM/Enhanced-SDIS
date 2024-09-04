function [pf, pf1, cov_t, n_cost, level_SuS,level_SDIS,cov_SuS,cov_SDIS] = SDIS_figure(g_fun,distr,nf,len,sigma,d,tarCoV)
%% Sequential directional importance sampling 
%{
Created by: Kai Cheng (kai.cheng@tum.de)
Based on: 
1. "Rare event estimation with sequential directional importance sampling", Cheng et al, Structural Safety, 100, 102291.
2. "Enhanced sequential directional importance sampling for structural reliability analysis", Cheng et al, In preparation.
3. "MCMC algorithms for subset simulation", Papaioannou et al, Probabilistic Engineering Mechanics 41 (2015) 83-103.
---------------------------------------------------------------------------
Input:
* nf     : important directions per level
* len    : length of each Markov chain
* g      : limit state function
* sigma  : initial simga
* d      : input dimension
* tarCoV : target coefficient of variation of the weights
---------------------------------------------------------------------------
Output:
* pf     : probability of failure
* pf1    : first probability of failure in SDIS
* cov_t  : coefficient of variation of SDIS estimator
* n_cost : total model evaluations
* level  : total number of intermediate levels
* level_SuS  : total number of intermediate levels in SuS for estimating pf1
* level_SDIS : total number of intermediate levels in SDIS for estimating pf/pf1
* cov_SuS    : coefficient of variation of SuS estimator for estimating pf1
* cov_SDIS   : coefficient of variation of SDIS estimator for estimating pf/pf1
%}

%% transform to the standard Gaussian space

if any(strcmp('Marginals',fieldnames(distr))) == 1   % use Nataf transform (dependence)
   n   = length(distr.Marginals);    % number of random variables (dimension)
   u2x = @(u) distr.U2X(u);          % from u to x
   
else   % use distribution information for the transformation (independence)
   % Here we are assuming that all the parameters have the same distribution !!!
   % Adjust accordingly otherwise
   n   = length(distr);                    % number of random variables (dimension)
   u2x = @(u)distr(1).icdf(normcdf(u));   % from u to x 

end

%% LSF in standard space

% g = @(u) g_fun(u2x(u)); 
 g = @(u) g_fun(u); 

%% Step 1 : Monte Carlo simulation 

% initilization
n                = 0;     % number of MCS population
n_failure        = 0;     % number of failure samples
model_evaluation = 0;     % number of model evaluations
sig              = sigma; % initial magnification factor of input standard deviation

mu = zeros(1,d); stdu = ones(1,d);         % mean and standard deviation of input variable

while n_failure < nf                       % sequential enrichment of initial MCS population

    n      = n + 1;   
    x(n,:) = lhsnorm(mu,diag(stdu.^2),1);  % generate random samples with MCS
    y(n)   = g(sig.*x(n,:));               % evaluate the auxiliary limit state function 

    if y(n) < 0
       n_failure  = n_failure + 1;
       ind(n_failure) = n;
    end
    
    if n > nf*10-1                         % 10% of the total samples
        break;
    end
    
end

if n_failure < nf
   b = prctile(y, 10);                     % intermediate failure threshold 
   ind = find(y < b);
else
   b = 0;
end
      
nf = length(ind);
model_evaluation = model_evaluation + n;  % update computational cost

%%  SuS
b
  figure
  subplot(1,2,1)
  g1 = @(x)normpdf(x);

  plot(x(:,1),x(:,2),'o','Markersize',5); hold on %Figure

  if d == 2
     xx    = -10:0.1:10;
     nnp   = length(xx);
     [X,Y] = meshgrid(xx);
     xnod  = cat(2,reshape(X',nnp^2,1),reshape(Y',nnp^2,1));

     Z     = g(3.*xnod); Z1 = g(xnod); Z2  = g(3.*xnod)-b(1);

     Z     = reshape(Z,nnp,nnp);
     Z1    = reshape(Z1,nnp,nnp);
     Z2    = reshape(Z2,nnp,nnp);

     Beta  =[1 2 3 4 5 6 7];
     v     = normcdf(-Beta);
     contour(Y,X,Z1,[0,0],'b-','LineWidth',1.5); hold on; % LSF
     contour(Y,X,Z,[0,0],'m-','LineWidth',1.5); hold on; % LSF
     contour(Y,X,Z2,[0,0],'g-','LineWidth',1.5); hold on; % LSF
     axis equal tight;
     axis([-10 10, -10 10])
  end

 % plot(X_original(:,1),X_original(:,2),'o','Markersize',5); hold on

 % plot(x_f(:,1),x_f(:,2),'o','Markersize',5); hold on %Figure

 xlabel('$u_1$',Interpreter="latex");
 ylabel('$u_2$',Interpreter="latex")
 title('First step',Interpreter="latex")
 % legend('MCS population','$G({\bf{u}})=0$','$G({3\bf{u}})=0$',Interpreter="latex")
 legend('MCS population','$G({\bf{u}})=0$','$G({3\bf{u}})=0$','$G({3\bf{u}})=41.26$',Interpreter="latex")

%% Sequential reduction of ---- SuS stage

if b > 0
   
   prob(1)     = length(ind)/n;                 % estimation of P_1 with MCS
   cov_sus(1)  = sqrt((1-prob(1))/n/prob(1));   % coefficient of variation of P_1       
 
   j = 1;   p0 = 0.1;   x_seeds = x(ind,:);  y_seeds = y(ind)';  nl = 10;   

   while b(j) > 0
    
       [x_mcmc, y_mcmc, ratio, alphak, indc] = MCMC_SuS(x_seeds,y_seeds, g,nf,10,sig,b(j)); % MCMC sampling
  
       model_evaluation = model_evaluation + length(y_mcmc)*9/10;    % update computational cost
  
       j = j + 1;
   
       if  indc < nf
           b(j) = 0;
       else
           b(j) = max(prctile(y_mcmc, p0*100),0);
       end
     
       % failure points in the next level

        ind = find(y_mcmc <= b(j)); 

       if b(j) > 0    
          
          x_seeds = x_mcmc(ind,:);
          y_seeds = y_mcmc(ind);
    
          ns  = length(y_mcmc);  % total sample size         
          nf1 = ns/nl;
 
          I_Fj     = reshape(y_mcmc <= b(j),nl,nf1);              % indicator function for the failure samples
          p_j      = (1/ns)*sum(I_Fj(:));                         % ~=p0, sample conditional probability
          gamma    = corr_factor(I_Fj,p_j,nl,nf1);                % corr factor 
          cov_sus(j) = sqrt(((1-p_j)/(ns*p_j))*(1+gamma));  

       else  % the last step

           y_remained = y_seeds(indc+1:end);
           x_remained = x_seeds(indc+1:end,:);
           
           ind_remained = find(y_remained < 0 );
           ny  = length(y_remained);      % seeds without move
           ns  = length(y_mcmc);          % total sample size         
           nf1 = ns/nl;
 
           I_Fj     = reshape(y_mcmc <= b(j),nl,nf1);                % indicator function for the failure samples
           p_j      = (sum(I_Fj(:)) + sum(y_remained<0))./(ns+ny);
           gamma    = corr_factor(I_Fj,p_j,nl,nf1);                  % corr factor (Ref. 2 Eq. 10)
          
           
           cov_sus(j) = sqrt((1-p_j)/p_j*(ns*(1+gamma) + ny))/(ns+ny);
           
           x_seeds = [x_mcmc(ind,:); x_remained(ind_remained,:)];
           y_seeds = [y_mcmc(ind);   y_remained(ind_remained)];

       end

       prob(j) = p_j;

   end

   pf1 = prod(prob);  cov_SuS = sqrt(sum(cov_sus.^2));

   x_f = x_seeds; y_f = y_seeds; 
   
else
    
   pf1     =  (nf-1)/(n-1);  
   cov_SuS = sqrt((1-pf1)/(n-2));
   x_f     = x(ind,:); 
   y_f     = y(ind); 
   prob    = pf1;

end

level_SuS = length(prob);

  subplot(1,2,2)
  g1 = @(x)normpdf(x);
 
  if d == 2
     xx    = -10:0.1:10;
     nnp   = length(xx);
     [X,Y] = meshgrid(xx);
     xnod  = cat(2,reshape(X',nnp^2,1),reshape(Y',nnp^2,1));

     Z     = g(3.*xnod); Z1 = g(xnod); 
     Z     = reshape(Z,nnp,nnp);
     Z1    = reshape(Z1,nnp,nnp);

     Beta  =[1 2 3 4 5 6 7];
     v     = normcdf(-Beta);
     contour(Y,X,Z1,[0,0],'b-','LineWidth',1.5); hold on; % LSF
     contour(Y,X,Z,[0,0],'m-','LineWidth',1.5); hold on; % LSF
     axis equal tight;
     axis([-10 10, -10 10])
  end

 % plot(X_original(:,1),X_original(:,2),'o','Markersize',5); hold on
 % plot(x(:,1),x(:,2),'o','Markersize',5); hold on %Figure
  plot(x_f(:,1),x_f(:,2),'o','Markersize',5); hold on %Figure

  xlabel('$u_1$',Interpreter="latex");
  ylabel('$u_2$',Interpreter="latex")
  for i=1:length(x_f)
     quiver(0,0,x_f(i,1),x_f(i,2), 'color', [0, 0.9, 0.9]);  hold on
  end
  title('Second step',Interpreter="latex")
  legend('G({\bf{u}})=0','G({3\bf{u}})=0','Failure samples','Importance directions',Interpreter="latex")


%% Sequential reduction of Sigma ---- SDIS stage
 
x0 = mu; y0 = g(x0);   % origin and its response

model_evaluation = model_evaluation + 1;    % update computational cost

 if b(1) > 0
    
   [y_u, ia, ic] = unique(y_f);
   
   [model_run, alpha, root, n_root] = Roots_Kriging(x_f(ia,:),y_f(ia),y0,g,sig); % find roots with Kriging model
   
    model_evaluation = model_evaluation + sum(model_run);  % update computational cost

    directions = alpha(ic,:);          % important direction     
    roots     = root(ic,:);
    n_roots   = n_root(ic);

 else
     
   [model_run, alpha, root, n_root] = Roots_Kriging(x_f,y_f,y0,g,sig);  % find roots with Kriging model

   model_evaluation = model_evaluation + sum(model_run);  % update computational cost

   directions = alpha;          % important direction     
   roots      = root;
   n_roots    = n_root;
     
 end
 
tarWk       = tarCoV;         % Target coefficient of variation of important weight
optimal_sig = sig;           
k = 1;

while optimal_sig > 1
  
   options     = optimoptions('fmincon','Display','off');            
   [optimal_sig, fval] = fmincon(@(x)Weight(x,sig(k),roots,tarWk,d,n_roots),1,[],[],[],[],1,sig(k),[],options);  % find optimal sigma
 
   if abs(optimal_sig-1) < 0.005
       optimal_sig = 1;
   end

   if fval > 0.01 && optimal_sig ~= 1 % refind optimal sigma

       ss = 1:(sig(k)-1)/100:sig(k);
       for i = 1:length(ss)
           [obj1(i), ~, ~ , ~] = Weight(ss(i),sig(k),roots,tarWk,d,n_roots);
       end

       [~,index] = min(obj1);
       [optimal_sig, ~] = fmincon(@(x)Weight(x,sig(k),roots,tarWk,d,n_roots),ss(index),[],[],[],[],1,sig(k),[],options);  % find optimal sigma
    
   end

   k = k + 1;  sig(k) = optimal_sig;

   [obj, weight, mu_w, std_w] = Weight(optimal_sig,sig(k-1),roots,tarWk,d,n_roots); % importance weight

   sk(k-1) = mu_w;                                     % mean of importance weight 
                                       
   cv(k-1) = sqrt((std_w/sk(k-1))^2/length(n_roots));   % coefficient of variation of importance weight  

   if optimal_sig == 1                     % convergence criterion
      break; 
   end
    
%%  Resampling   

  [seeds, ind]= Resample(weight,roots, n_roots, directions, sig, k, d, nf);      % resampling 
  
  [x_mcmc, y_mcmc, ratio, ind1] = MCMC_SDIS(seeds, g, nf, len, sig(k));            % MCMC for sampling

  % [x_mcmc, y_mcmc, ind1] = ESS_SDIS(seeds, g, nf, len, sig(k));            % MCMC for sampling

  ind2 = 1:1:nf;  ind2(ind1) = [];                                               % index of markov chains which all samples are rejected
     
  x_stable = x_mcmc(len:len:end,:);  y_stable = y_mcmc(len:len:end);             % select the desired stable samples
    
  model_evaluation = model_evaluation + nf*len;                                  % update computational cost

  [model_run, direction1, root1, n_root1] = Roots_Kriging(x_stable(ind1,:),y_stable(ind1),y0,g,sig(k));  % find the roots of the failure samples
   
   direction2 = zeros(nf,d);
   root2      = zeros(nf,20);
   n_root2    = zeros(1,nf);

   % new directions and roots
   direction2(ind1,:) = direction1;  
   root2(ind1,:)      = root1;  
   n_root2(ind1)      = n_root1; 

   % old directions and roots
   direction2(ind2,:) = directions(ind(ind2),:); 
   root2(ind2,:)      = roots(ind(ind2),:)*sig(k-1)/sig(k);
   n_root2(ind2)      = n_roots(ind(ind2)); 
   
   directions = direction2;
   roots      = root2;
   n_roots    = n_root2;
   
   model_evaluation = model_evaluation + sum(model_run);                                % update computational cost

 end

 pf       = pf1*prod(sk);                     % failure probability 
 cov_SDIS = sqrt(sum(cv.^2));                 % coefficient of variation of SDIS
 cov_t    = sqrt(cov_SuS.^2 + cov_SDIS.^2);   % total coefficient of variation 
 n_cost   = model_evaluation;                 % total model evaluations
 level_SDIS  = length(sk);                    % total number of intermediate levels

end




