clear
close all;

%%% Multi objective Kalman Filter for a simple example
%CHECK the interval a,b and the initial datum u

%% Initial data
a=-2;
b=4;      %domain

Nsample=10;  %number of the ensembles
K=10;    % number of the observations (sensors :))
NL=20;         %Lambda discretization parameter
d=1;          %dimension of "u"
Nmax_it=100;
lambda=linspace(0.0,1,NL);

%u = a + (b-a).*rand(d,d);
%u=1;
mu=zeros(K,1);
Gamma=eye(K);%0.05*eye(K);
noise = mvnrnd(mu, Gamma);

y1= zeros(K,1);%G1_scalar(u, K) + noise';    
y2= zeros(K,1);%G2_scalar(u, K) + noise';   % should be a K x 1 vector 


Gamma_inv=Gamma \ eye(K);
%plot (y1,'*'); hold on; plot (y2,'+')

%% EnKF procedure
rng(2);  %%%%%% BE CAREFULLLL!!!!!!!!!!!!!!!!
u0 = a + (b-a)*(rand(1,Nsample));    %u+noise';%a + (b-a).*(rand(Nsample,1));
%u0 = sort(u0);
plot(u0,'*');
figure
d = 1;
u_final=zeros(d, Nsample, NL);
us = zeros(d, NL);
N = zeros(1, NL);
xax = zeros(1, NL);
yax = zeros(1, NL);
g = zeros(K, NL);
g1 = zeros(K, NL);
g2 = zeros(K, NL);
df1 = zeros (K,NL);
df2 = zeros (K,NL);

% Initial data for moment equations

m0=sum(u0,2)/Nsample;               % Initial momentum
E0=sum(u0.*u0,2)/Nsample;             % initial energy
ml0=0;
El0=0;

um=u0;

for i=1:NL
 %  [us(i), unp1(:,i), uit(:,i)]=EnKF(lambda(i),y1,y2,Nsample,u0,Gamma,Gamma_inv,M);
   [us(:, i), u_final(:,:,i), N(i)]=EnKF(lambda(i),y1,y2,Nsample,um,Gamma,Gamma_inv,K,Nmax_it);

   g(:, i)=G_scalar(lambda(i),us(:,i), K);
   g1(:, i)=G1_scalar(us(:,i), K);
   g2(:, i)=G2_scalar(us(:,i), K);

   % Pareto stuff
   df1(:,i) = g1(:,i) - y1(:);%y(y1,y2,lambda(i)) - G1_scalar(us(:, i), K);
   xax(i) = df1(1,i);%0.5 * dy' * Gamma_inv * dy;
   
   df2(:,i) = g2(:,i) - y2(:);%y(y1,y2,lambda(i)) - G2_scalar(us(:, i), K);
   yax(i) = df2(1,i);%0.5 * dy' * Gamma_inv * dy; 
   
% preparation for the next step

t_interval = [0 10];
init_cond = [m0,0,E0,0]';

[t,mom] = ode45(@(t,Y) ode_system(t,Y,lambda(i),K,y1,y2), t_interval , init_cond);

%new initial data
u_eps=rand(1,Nsample);
%um=u0;
um=mom(size(mom(:,1),1),1)+u_eps;
end


%us(1)-u
%plot(uit(:)-us(1))
%% Plot
 plot(us,g1(1,:),'*-b'); hold on;
 plot(us,g2(1,:),'o-r');
 plot(us,g(1,:),'-.k');
 legend('g1','g2','g');
%  figure
% d=sqrt(xax.^2+yax.^2);
% [r, idx] = min(d);
% phi = linspace(0,2.1*pi, 100);
% xcirc = r*cos(phi);
% ycirc=r*sin(phi);
figure 
plot(t,mom(:,1),'b',t,mom(:,2),'r',t,mom(:,3),'k',t,mom(:,4),'g','LineWidth',2); legend('m','m_\lambda','E','E_lambda')

figure
title('Pareto')
for i=1:NL
    if (N(i)<Nmax_it)
        plot (xax(i),yax(i),'o-r','LineWidth',1);hold on;
    else
        plot (xax(i),yax(i),'o-k','LineWidth',1);hold on;
    end
end
title('Pareto')
xlabel('y_1-G_1')
ylabel('y_2-G_2')
hold on
exact_Pareto;
%plot(xcirc, ycirc, 'r')
%save('Pareto_point.mat','xax','yax')
save('Naive_point.mat','xax','yax')

figure
plot (lambda , xax);
xlabel('lambda')
%ylabel('y_1-G_1)')
hold on;

%figure
plot (lambda, yax);
xlabel('lambda')
%ylabel('norm(y_2-G_2)')
