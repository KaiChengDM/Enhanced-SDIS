
function [seeds, ind] = Resample(weight,roots, n_roots, direction, sig, k, d, nchain)

%% Resample for important directions

 nsamlev     = length(weight);                              % number of samples
 ind         = randsample(nsamlev,nchain,true,weight);      % resampling of directions 
 
 % weight    = weight(ind);
 root_seeds  = roots(ind,:);                                % resampled roots along important directions
 n_seeds     = n_roots(ind); 
 
%% Resample of radius

 for i = 1 : nchain   
    
     r_accepted = [];

     root = sig(k-1)./sig(k)*root_seeds(i,:);

     while isempty(r_accepted)

           r   = sqrt(chi2rnd(d,1,10^3));    % random radius 

           if n_seeds(i) > 1

              if mod(n_seeds(i),2)  == 0     % even case

                 for j = 1:2:n_seeds(i)    

                    r_accepted = [r_accepted r(root(j+1) > r & r > root(j))] ; 

                 end

              else                            % odd case

                  for j = 1:2:n_seeds(i)-1 

                      r_accepted = [r_accepted r(root(j+1) > r & r > root(j))] ; 

                  end

                  r_accepted = [r_accepted r(r > root(n_seeds(i)))] ;

              end

           else
               r_accepted = [r_accepted r(r > root(1))];  % accepted radius 
           end

     end

    % r_accepted = resample_radius(root_seeds(i,:),n_seeds(i),sig,k,d);
    
    num  = size(r_accepted,2);     % number of accepted radius
   
    ind1 = randsample(num,1);      % resampling of radius

    seeds(i,:) = r_accepted(ind1).*direction(ind(i),:);
    
    % seeds(i,:) = r_accepted.*direction(ind(i),:);  % transform radius and direction into Cartesian coordinates 

 end

end



