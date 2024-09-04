function [roots, m_run] = Kriging(g1,inputpar,d,n_size,sigma)
 
 % addpath(genpath('Kriging-and-GE-Kriging-model-toolbox-main'));
 
 addpath(genpath('dace'));

 alpha = 10^-10;  lb = sqrt(chi2inv(alpha/2,d))/sigma; ub = sqrt(chi2inv(1-alpha/2,d));  % confidence interval

 x = inputpar.x;  y = inputpar.y; 

 len = ub - lb;  

 if abs(x(2)-ub) < len/3
     x(3) = (lb+x(2))/2; y(3) = g1(x(3));
     if n_size > 3                 % 4 initial training samples 
         if y(3) > 0
            x(4) = (x(3)+x(2))/2;   
         else
            x(4) = (lb+x(3))/2;   
         end
         y(4) = g1(x(4));
     end
 else
    x(3) = (ub+x(2))/2;

    y(3) = g1(x(3));

     if n_size > 3                  % 4 initial training samples 
        x(4) = (lb+x(2))/2; y(4) = g1(x(4));    
     end
     
 end
 
 options = optimoptions('fsolve','Display','off','OptimalityTolerance',10^-10,'FunctionTolerance',10^-10,'MaxIterations',400);
 
 options1 = optimoptions('fmincon','Display','off');

 num = 20;  seeds = (lb:(ub-lb)/(num-1):ub)';
 
 d     = 1;
 theta = 0.1*ones(1,d); 
 lob   = 10^-4*ones(1,d);
 upb   = 20*ones(1,d);

while 1

    try
       [Kriging_model , perf] = dacefit(x,y,@regpoly1,@corrmatern,theta,lob,upb); 
    catch
       [Kriging_model , perf] = dacefit(x,y,@regpoly0,@corrmatern,theta,lob,upb);
    end

    F = @(t)predictor(t,Kriging_model);

    root = [];
   
    for i = 1:num   
        [root(i), fval(i), exitflag(i)] = fsolve(F,seeds(i),options);  % find root  
    end

    x0 = (lb:(ub-lb)./(20-1):ub)';       % with 10 initial points

    for i = 1:20
        [x_best(i), EDR_value(i)]= fmincon(@(r)EDR_function(r,F),x0(i),[],[],[],[],lb,ub,[],options1);
    end

     [EDR_max, ind] = min(EDR_value);
     
     if  abs(EDR_max)/mean(abs(y)) < 0.0005 % convergence criterion
         break;
     end

     x_next = x_best(ind);
   
     x = [x;  x_next];  y = [y; g1(x_next)];
     
%       r = (0:(ub+2-0)./(500-1):ub+2)';
% 
%  for i = 1:length(r)
%      yt(i) = g1(r(i));
%      mu(i) = F(r(i));
%  end
% 
% 
% if n_size < 4
%    figure
%    plot(r,yt','LineWidth',2); hold on
%    plot(r,mu,'LineWidth',2); hold on
%    plot(x(1:3),y(1:3),'*','LineWidth',1); hold on
%    plot(lb,g1(lb),'square','LineWidth',1); hold on
%    plot(ub,g1(ub),'square','LineWidth',1); hold on
%    plot(x(4:end),y(4:end),'o','LineWidth',1); hold on
%    plot(r,zeros(1,length(r)),'g--','LineWidth',1); hold on
%    legend('True limit state','Kriging mean','Initial samples','lb','ub','Enriched samples')
% else
%    figure
%    plot(r,yt','LineWidth',2); hold on
%    plot(r,mu,'LineWidth',2); hold on
%    plot(x(1:4),y(1:4),'*','LineWidth',1); hold on
%    plot(lb,g1(lb),'square','LineWidth',1); hold on
%    plot(ub,g1(ub),'square','LineWidth',1); hold on
%    plot(x(5:end),y(5:end),'o','LineWidth',1); hold on
%    plot(r,zeros(1,length(r)),'g--','LineWidth',1); hold on
%    legend('True limit state','Kriging mean','Initial samples','lb','ub','Enriched samples')
% end


end
   
   root = root(abs(F(root')) < 0.01);       % double check fake roots
     
   root = unique(round(root,4))';           % remove repeated roots

   ind  = find((lb < root) & (root < ub));  % find roots within the confidence interval
   
   if length(ind) < 1                       % no roots within the confidence interval    
       if F(lb) > 0 
          roots = ub;
       else
          roots = lb;
       end
   else
       if F(lb) > 0 
          roots = root(ind);
       else
          roots = [lb; root(ind)];
       end
   end
   
   roots = sort(roots,'ascend');
   
   m_run = length(y) - 2; % model evaluations on a direction
   
 
  %% plot
 
%  r = (0:(ub+2-0)./(500-1):ub+2)';
% 
%  for i = 1:length(r)
%      yt(i) = g1(r(i));
%      mu(i) = F(r(i));
%  end
% 
% 
% if n_size < 4
%    figure
%    plot(r,yt','LineWidth',2); hold on
%    plot(r,mu,'LineWidth',2); hold on
%    plot(x(1:3),y(1:3),'*','LineWidth',1); hold on
%    plot(lb,g1(lb),'square','LineWidth',1); hold on
%    plot(ub,g1(ub),'square','LineWidth',1); hold on
%    plot(x(4:end),y(4:end),'o','LineWidth',1); hold on
%    plot(r,zeros(1,length(r)),'g--','LineWidth',1); hold on
%    legend('True limit state','Kriging mean','Initial samples','lb','ub','Enriched samples')
% else
%    figure
%    plot(r,yt','LineWidth',2); hold on
%    plot(r,mu,'LineWidth',2); hold on
%    plot(x(1:4),y(1:4),'*','LineWidth',1); hold on
%    plot(lb,g1(lb),'square','LineWidth',1); hold on
%    plot(ub,g1(ub),'square','LineWidth',1); hold on
%    plot(x(5:end),y(5:end),'o','LineWidth',1); hold on
%    plot(r,zeros(1,length(r)),'g--','LineWidth',1); hold on
%    legend('True limit state','Kriging mean','Initial samples','lb','ub','Enriched samples')
% end


end

