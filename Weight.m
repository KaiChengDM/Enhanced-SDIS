function [obj, weight, mu_w, sig_w] = Weight(x,sig,root,tarWk,d,n_roots)

 % Compute the variability of the important weight

  weight = pf_directional(root*sig/x,d,n_roots)./pf_directional(root,d,n_roots);

  mu_w = mean(weight); sig_w = std(weight);

  obj = abs(sig_w/mu_w- tarWk);   
        
  
end

