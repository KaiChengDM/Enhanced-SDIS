
clc; clear; 
%% preparation

d = 2;          % number of dimensions

kappa = 1;

Beta1  = 3.5;

g_fun{1} = @(x) -sum(x')'/sqrt(d) + Beta1 ;

g_fun{2} = @(x)-sum(x')'/sqrt(d) + Beta1 - kappa*(x(:,1)-x(:,2)).^2/10;    % limit state function   %g  = @(x)-sum(x')'/sqrt(d) + Beta;

g_fun{3} = @(x)-sum(x')'/sqrt(d) + Beta1 + kappa*(x(:,1)-x(:,2)).^2/10;    % limit state function  

%% Figure 1
% 
mu = zeros(1,d); stdu = eye(d);         % mean and standard deviation of input variable
n_mcs = 1000;

g1 = @(x)normpdf(x);

for i = 1:3

g = g_fun{i};

subplot(1,3,i)

Z = []; Z1 = [];  Z2 = []; 

if d == 2

     xx    = -10:0.1:10;
     nnp   = length(xx);
     [X,Y] = meshgrid(xx);
     xnod  = cat(2,reshape(X',nnp^2,1),reshape(Y',nnp^2,1));

     Z  = g(xnod); Z1 = g(2*xnod); Z2 = g(3*xnod);  

     Z3 = g1(xnod);  

     Z     = reshape(Z,nnp,nnp);
     Z1    = reshape(Z1,nnp,nnp);
     Z2    = reshape(Z2,nnp,nnp);
     Z3    = reshape(prod(Z3')',nnp,nnp);

     u = lhsnorm(mu,stdu,n_mcs);
     plot(u(:,1),u(:,2),'bo');hold on;

     Beta  =[1 2 3 4 5 6 7];
     v     = normcdf(-Beta);

     contour(Y,X,Z,[0,0],'b-','LineWidth',1.5); hold on; % LSF
     contour(Y,X,Z1,[0,0],'r-','LineWidth',1.5); hold on; % LSF
     contour(Y,X,Z2,[0,0],'m-','LineWidth',1.5); hold on; % LSF
     % contour(Y,X,Z2,[b,b],'y-','LineWidth',1.5); hold on; % LSF
     % contour(Y,X,Z3,v,'g:','LineWidth',1); hold on; % LSF
     % contour(Y,X,Z3,5,'LineWidth',1); hold on; % LSF

     axis equal tight;
     axis([-10 10, -10 10])
  end

  legend('Samples','$G({\bf{u}})=0$','$G(2{\bf{u}})=0$','$G(3{\bf{u}})=0$','Contour of input PDF',Interpreter="latex")

  xlabel('$u_1$',Interpreter="latex");
  ylabel('$u_2$',Interpreter="latex")
  hold on
  title('Linear LSF',Interpreter="latex")

end
% 
%% Figure 2

d = 2;

% g = @(x) -sum(x')'/sqrt(d) + Beta1  + kappa*(x(:,1)-x(:,2)).^2/0.3;

g  = @(x)30./((4*(x(:,1)+2).^2/9+x(:,2).^2/25).^2 + 1) + 20./(((x(:,1)-2.5).^2/4 + (x(:,2)-0.5).^2/25).^2+1) - 5;

mu = zeros(1,d); stdu = eye(d);  n_mcs = 1000;

u = lhsnorm(mu,stdu,n_mcs);

y = g(3.*u);

b = prctile(y, 10);

g1 = @(x)normpdf(x);

figure

Z = []; Z1 = [];  Z2 = []; 

if d == 2
     xx    = -10:0.1:10;

     nnp   = length(xx);

     [X,Y] = meshgrid(xx);

     xnod  = cat(2,reshape(X',nnp^2,1),reshape(Y',nnp^2,1));

     Z  = g(xnod); Z1 = g(3*xnod); Z2 = g(3*xnod);  

     Z3 = g1(xnod);  

     Z     = reshape(Z,nnp,nnp);

     Z1    = reshape(Z1,nnp,nnp);
     Z2    = reshape(Z2,nnp,nnp);
     Z3    = reshape(prod(Z3')',nnp,nnp);

     u = lhsnorm(mu,stdu,n_mcs);

     plot(u(:,1),u(:,2),'bo');hold on;

     Beta  =[1 2 3 4 5 6 7];
     v     = normcdf(-Beta);

     contour(Y,X,Z,[0,0],'b-','LineWidth',1.5); hold on; % LSF
     contour(Y,X,Z2,[0,0],'m-','LineWidth',1.5); hold on; % LSF
     contour(Y,X,Z1,[b,b],'g-','LineWidth',1.5); hold on; % LSF
     % contour(Y,X,Z3,v,'g:','LineWidth',1); hold on; % LSF

     axis equal tight;
     axis([-10 10, -10 10])

  end

 legend('Samples','$G({\bf{u}})=0$','$G(3{\bf{u}})=0$','$G(3{\bf{u}})=4.28$','Contour of input PDF',Interpreter="latex")

 xlabel('$u_1$',Interpreter="latex");
 ylabel('$u_2$',Interpreter="latex")
 hold on
 title('Linear LSF',Interpreter="latex")