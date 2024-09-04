
function  [x_mcmc,y_mcmc] = ESS_SDIS(seeds,g, nf,len,sigma)

%% MCMC sampling: Elliptical slice algorithm

  [m, n]      = size(seeds);   % size of seeds
  seedslength = nf;            % number of Markov chains 
  lenchain    = len;           % length of each Markov chain
  alphak = zeros(m,1);         % space for the standard deviation

  ii = 0; ind = [];            % num = 0;

  for k = 1 : seedslength
  
      x_current = seeds(k,:);   % Initial state
      y_current = 0;            % dummy value 0
      flag = 0 ;
      
      theta = rand*2*pi; 

      lb = theta - 2*pi;
      ub = theta;

      u = mvnrnd(zeros(n,1), eye(n),1);

      while 1

         x_candidate = x_current.*cos(theta) + u.*sin(theta);
         
         y_candidate = g(sigma.*x_candidate);    
     
         if y_candidate < 0

            ii = ii + 1;
            x_current = x_candidate; 
            y_current = y_candidate;

            x_mcmc(ii,:) = x_current;
            y_mcmc(ii,:) = y_current;

            break;

         else

            % lb =  lb*exp(-1);
            % ub =  ub*exp(-1);

            if theta < 0
                lb = theta;
            else
                ub = theta;
            end

            theta = lb+rand*(ub-lb);

         end
 
     end
   
  end
    
end


