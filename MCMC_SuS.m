
function  [x_mcmc, y_mcmc, alpha, alphak,k] = MCMC_SuS(x_seeds,y_seeds, g,nf,nl,sigma,b)

%% MCMC sampling: conditional sampling M-H (CSM-H) algorithm

  [m, n]     = size(x_seeds);  % size of seeds
  seedsnum    = m;             % number of Markov chains 
  lenchain    = nl;            % length of each Markov chain
  beta        = 0.5;           % initial beta    
  alphak      = zeros(m,1);    % space for the standard deviation
  sigmaf      = 1;             % initial standard deviation
  lambda      = 0.6;           % initial lambda
  counta      = 0; 
 
  ii = 0; 

 for k = 1 : seedsnum

     x_current = x_seeds(k,:);   % initial state
     y_current = y_seeds(k);     % initial state

     ii = ii + 1;
     
     x_mcmc(ii,:) = x_current;
     y_mcmc(ii,:) = y_current;

     for j = 1:nl-1
       
         rhok = beta;
         x_candidate = normrnd(rhok*x_current', sqrt(1-rhok^2))';   % candidate state
         y_candidate = g(sigma.*x_candidate);    
  
         if  y_candidate < b                 % accept or reject
             x_current = x_candidate;
             y_current = y_candidate;
             alphak(k) = alphak(k)+1/lenchain;
         end

         ii = ii + 1;
         x_mcmc(ii,:) = x_current;
         y_mcmc(ii,:) = y_current;

     end

      nfy = sum(y_seeds(k+1:end)<0);

      nf_mcmc = sum(unique(y_mcmc)<0);  % for the last step
   
      if nf_mcmc > nf-nfy
          break; 
      end

     adapchains = 10;
     % check whether to adapt now
     if mod(k,adapchains) == 0 
        % num = num +1;
        % mean acceptance rate of last adap_chains
        alpha_mu = mean(alphak(k-adapchains+1:k));
        counta   = counta + 1;
        gamma    = counta^(-0.5);
        lambda   = exp(log(lambda)+gamma*(alpha_mu-0.44));    
        % compute parameter rho
        sigmafk = min(lambda*sigmaf,1);
        beta   = sqrt(1-sigmafk.^2);
     end
     
 end

  alpha = mean(alphak(1:k));
  
end


