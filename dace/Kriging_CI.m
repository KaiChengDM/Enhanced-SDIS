clc; clear all; 

% g = @(x)3.*(1-x(:,1)).^2.*exp(-x(:,1).^2-(x(:,2)+1).^2)-10.*(0.2.*x(:,1)-x(:,1).^3-x(:,2).^5).*exp(-x(:,1).^2-x(:,2).^2)-1./3.*exp(-(x(:,1)+1).^2-x(:,2).^2);

g = @(x) x(:,1).*sin(x(:,1));

d = 1; N = 7;    % n-Initial samples size; n1-test samples size;

Lb = 0; Ub = 15;

pp = sobolset(d,'Skip',3); u=net(pp,N);  

for i = 1:d
  x(:,i)=u(:,i)*(Ub(i)-Lb(i))+Lb(i);
end

y=g(x); 

%% Kriging
addpath(genpath('dace'));

theta =1.*ones(1,d);lob=10^-3.*ones(1,d);upb=10.*ones(1,d);%%%%theta??????????????????

[dmodel, perf]=dacefit(x,y,@regpoly0,@corrgauss,theta,lob,upb);%%????????????

F = @(t)predictor(t,dmodel);

xx = Lb :0.01 :Ub;

[yy, vv]= F(xx'); yy1 = g(xx');

up = yy + 1.96.*sqrt(vv);
lp = yy - 1.96.*sqrt(vv);

figure
plot (x,y,'ro','LineWidth',2); hold on
plot (x(7),y(7),'m*','LineWidth',2.5); hold on

plot (xx,yy1,'b-','LineWidth',1.5); hold on
plot (xx,yy,'b--','LineWidth',1.5);  hold on
plot (xx,up,'r:','LineWidth',1.5);  hold on
plot (xx,lp,'r:','LineWidth',1.5); hold on


for i = 1 :length(lp)
   yy2 = lp(i):0.01: up(i);
   xx2 = xx(i).*ones(length(yy2),1);
   plot (xx2,yy2,'g','LineWidth',1.5); hold on
end

plot (x,y,'ro','LineWidth',2); hold on
plot (x(7),y(7),'mo','LineWidth',2.5); hold on

plot (xx,yy1,'b-','LineWidth',1.5); hold on
plot (xx,yy,'b--','LineWidth',1.5);  hold on
plot (xx,up,'r:','LineWidth',1.5);  hold on
plot (xx,lp,'r:','LineWidth',1.5); hold on

xlabel('x','Fontsize',15);
ylabel('y','Fontsize',15)
 legend('Samples','Added sample','True response','Kriging predictor',' Conf. interval')
%legend('Samples','True response','Kriging predictor',' Conf. interval')

 [value, location] = max(vv);

 x = [x ; xx(location)];
 y = [y ; g(xx(location))];

