clear
close all;

%%% Multi objective Kalman Filter for a simple example
%CHECK the interval a,b and the initial datum u
%Select 1 for direct approach 
%select 2 for adaptive approach

%strategy=1; %Direct approach 

strategy=2;  %Adaptive approach


%% Initial data
% domain
a=-1;
b=1;     

Nsample=20;       %number of the ensembles
K=10;             % number of the observations 
NL=34;            %Lambda discretization parameter
d=1;              %dimension of "u"
Nmax_it=1000;
Nmax_it_lambda=400;
cost_up=0.0001;

mu=zeros(K,1);
Gamma=eye(K);
noise = mvnrnd(mu, 0.001*Gamma);

y1= zeros(K,1)+ noise'; 
y2= zeros(K,1)+noise';


Gamma_inv=Gamma \ eye(K);

%% EnKF procedure

rng(2);  % For reproducibility 
u0 = a + (b-a)*(rand(1,Nsample));    
plot(u0,'*');       %plot of the initial ensemble

%initialization of the vectors
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

i=1;
lambda(i)=0.0;

if strategy==1   %Direct approach
    lambda=linspace(0.0,1.005,NL);
    for i=1:NL
        [us(:, i), u_final(:,:,i), N(i),phi(i)]=EnKF(lambda(i),y1,y2,Nsample,um,Gamma,Gamma_inv,K,Nmax_it);

         g(:, i)=G_scalar(lambda(i),us(:,i), K);
         g1(:, i)=G1_scalar(us(:,i), K);
         g2(:, i)=G2_scalar(us(:,i), K);
 
         % Pareto front
         df1(:,i) = g1(:,i) - y1(:);
         xax(i) = df1(1,i);

         df2(:,i) = g2(:,i) - y2(:);
         yax(i) = df2(1,i);    
    end
    save('Naive_point.mat','xax','yax')
    NLA=i+1;

else

    while(lambda(i)<1-0.1 && i<Nmax_it_lambda) %Adaptive appoach

       [us(:, i), u_final(:,:,i), N(i),phi(i)]=EnKF(lambda(i),y1,y2,Nsample,um,Gamma,Gamma_inv,K,Nmax_it);

       g(:, i)=G_scalar(lambda(i),us(:,i), K);
       g1(:, i)=G1_scalar(us(:,i), K);
       g2(:, i)=G2_scalar(us(:,i), K);

       % Pareto front
       df1(:,i) = g1(:,i) - y1(:);
       xax(i) = df1(1,i);

       df2(:,i) = g2(:,i) - y2(:);
       yax(i) = df2(1,i);

       % preparation for the next step

        t_interval = [0 8];
        m0=sum(um,2)/Nsample;     %Initial momentum
        init_cond = [m0,0,E0,0]';
  
        [t,mom] = ode45(@(t,Y) ode_sys(t,Y,lambda(i),K,y1,y2,m0), t_interval , init_cond);

    %new initial data 
    s = rng;
    um = normrnd(mom(size(mom(:,1),1),1),mom(size(mom(:,1),1),3),[1,Nsample]);
    moment_in=mom(size(mom(:,1),1),1);
    lambda(i+1)=lambda(i)+cost_up/abs(mom(size(mom(:,1),1),2));

    if (mod(i,50)==0) %For checking
         i
     end

    i=i+1;

    end

save('Pareto_point.mat','xax','yax')

end
NLA=i-1  %number of lambda in the adaptive procedure


%% Plot
 fs=15;
 lw=1.5;
 figure
 title('Pareto')
 set(gca,'FontSize',fs);
 hold on
 %exact_Pareto;
 load('ParetoFronts.mat');  %Data for exact Pareto front for the toy examples
 scatter(Pf1(:,1),Pf1(:,2),'r.'); 
 for i=1:NLA
    if (N(i)<Nmax_it)
        plot (xax(i),yax(i),'* r','LineWidth',2);hold on;
    else
        plot (xax(i),yax(i),'* k','LineWidth',2);hold on;
    end
end
title('Pareto')
xlabel('y_1-G_1')
ylabel('y_2-G_2')
set(gca,'FontSize',fs);

figure
histogram(lambda(1:NLA+1),NLA+1)
set(gca,'FontSize',15);
title('Distribution of \lambda')
xlabel('Values of \lambda')
ylabel('Frequency')
axis([-0.05 1.05 0 1.2])
