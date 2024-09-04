function pf = pf_directional(root,d,n_roots)

 % Compute the failure probability on every important direction
 
 nf = size(root,1); 
 
 for i = 1 : nf
     
      n_root = n_roots(i);
      
      if n_root == 1                                                                 % one root 
         pf(i) = (1-chi2cdf(root(i,1).^2,d)); 
      elseif mod(n_root,2) == 0                                                      % even
         pf(i) = ((-1).^(1:n_root)*chi2cdf(root(i,1:n_root).^2,d)');  
      else                                                                           % odd
         pf(i) = (1-chi2cdf(root(i,n_root).^2,d) + (-1).^(1:n_root-1)*chi2cdf(root(i,1:n_root-1).^2,d)');
      end

 end

end

