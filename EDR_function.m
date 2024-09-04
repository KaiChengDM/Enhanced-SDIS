function EDR = EDR_function(r,F)

% Compute the Expected risk function

  % [mu, mse] = F(r);
 
 [mu, ~ ,mse] = F(r);

 s = sqrt(real(mse));

 EDR = -(-sign(mu).*mu.*normcdf(-sign(mu).*mu./s,0,1) + s.*normpdf(mu./s,0,1));
  
end

