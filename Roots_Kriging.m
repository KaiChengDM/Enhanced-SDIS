function [model_run, a, root, n_root] = Roots_Kriging(x_failure,y_failure,y0,g,sigma)

%% Find the root on specific direction

[m, n] = size(x_failure);     % failure sample size
root   = zeros(m,20);         % store the roots
n_size = 3;
n_root = 1;

for i = 1 : m

   e      = diag(ones(1,n));
   norm_x = norm(x_failure(i,:));
   a(i,:) = x_failure(i,:)./norm_x;  % direction defined by failure samples
   
   e(:,1) = a(i,:)';    
   [v, R] = qr(e);                  % define the rotatation matrix
   v(:,1) = a(i,:)';
  
   g0 = @(x)g(sigma*x*v');          % rotated performance function
   
   g1 = @(u)g0([u zeros(1,n-1)]);   % one-dimensional performance function
 
   inputpar.x = [0; norm_x];
   
   inputpar.y = [y0; y_failure(i)];

   [root_i, m_run(i)] = Kriging(g1,inputpar,n,n_size,sigma);  % find root 
   
   if mean(m_run) > 2  % reset the initial sample size
       n_size = 4;
   else
       n_size = 3;
   end

    n_root(i) = length(root_i); 
  
    root(i,1:n_root(i)) = root_i; 
  
    model_run(i) = m_run(i);  % computational cost 

end

end

