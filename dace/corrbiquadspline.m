function  [r, dr1 ]= corrbiquadspline(theta, d)
%CORRSPLINE  Biquadratic spline correlation function,
%
%           n
%   r_i = prod S(theta_j*d_ij) ,  i = 1,...,m
%          j=1
%
% with
%           1 - 15x^2 + 30x^3   for   0 <= x <= 0.5
%   S(x) =  1.25(1 - x)^3       for  0.5 < x < 1
%           0                   for    x >= 1
% If length(theta) = 1, then the model is isotropic:
% all  theta_j = theta.
%
% Call:    r = corrspline(theta, d)
%          [r, dr] = corrspline(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation
% dr    :  m*n matrix with the Jacobian of r at x. It is
%          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
%          where S(i,:) is the i'th design site. 

% hbn@imm.dtu.dk  
% Last update May 30, 2002

[m n] = size(d);  % number of differences and dimension of data
if  length(theta) == 1
  theta = repmat(theta,1,n);
elseif  length(theta) ~= n
  error(sprintf('Length of theta must be 1 or %d',n))
else
  theta = theta(:).';
end
mn = m*n;   ss = zeros(mn,1);
xi = reshape(abs(d) .* repmat(theta,m,1), mn,1);

% Contributions to first and second part of spline
 
i1 = find(xi <= 0.4);
i2 = find(0.4 < xi & xi < 1);

if  ~isempty(i1)
  ss(i1) = 1 - 15*xi(i1).^2 + 35*xi(i1).^3 - 195/8*xi(i1).^4;
end
if  ~isempty(i2)
  ss(i2) = 5/3 - 20/3*xi(i2) + 10*xi(i2).^2 - 20/3*xi(i2).^3 + 5/3*xi(i2).^4;
end

ss = reshape(ss,m,n);

r = prod(ss, 2);

%%
if  nargout > 1  % get Jacobian

  u = reshape(sign(d) .* repmat(theta,m,1), mn,1);

  dr = zeros(mn,1);
  if  ~isempty(i1)
    dr(i1) = -u(i1) .* ( -30*xi(i1) + 105*xi(i1).^2 - 195/2*xi(i1).^3);
  end
  if  ~isempty(i2)
    dr(i2) = -u(i2) .* (-20/3 + 20*xi(i2)- 20*xi(i2).^2 + 20/3*xi(i2).^3);
  end

  ii = 1 : m; dr1=dr;
  for  j = 1 : n
    sj = ss(:,j);  ss(:,j) = dr1(ii);
    dr1(ii) = prod(ss,2);
    ss(:,j) = sj;   ii = ii + m;
  end

  dr1 = reshape(dr1,m,n);

end

end  % Jacobian