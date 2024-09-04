
function  [x_mcmc,y_mcmc, alpha, ind ] = MCMC_SDIS(seeds,g, nf,len,sigma)

%% MCMC sampling: conditional sampling M-H (CSM-H) algorithm

  [m, n]      = size(seeds);   % size of seeds
  seedslength = nf;            % number of Markov chains 
  lenchain    = len;           % length of each Markov chain
  beta   = 0.7;                % initial beta    
  alphak = zeros(m,1);         % space for the standard deviation
  sigmaf = 1;                  % initial standard deviation
  counta = 0; 
  lambda = 0.6;
 
  ii = 0; ind = [];            % num = 0;

  for k = 1 : seedslength
  
      x_current = seeds(k,:);   % Initial state
      y_current = 0;            % dummy value 0
      flag = 0 ;
      
      for i = 1 : lenchain

         rhok = beta;
         x_candidate = normrnd(rhok*x_current', sqrt(1-rhok^2))';
         y_candidate = g(sigma.*x_candidate);    
   
         ii = ii + 1;
  
         if y_candidate < 0
            x_current = x_candidate; 
            y_current = y_candidate;
            alphak(k) = alphak(k)+1/lenchain;
            flag = 1;
         end

         x_mcmc(ii,:) = x_current;
         y_mcmc(ii,:) = y_current;
      
     end
   
     if flag == 1        % at least one sample is accepted
        ind = [ind, k];  % index of Markov chain that all samples are rejected
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
  
  alpha = mean(alphak);
  
end


