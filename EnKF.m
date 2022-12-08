function [mea,unp1,N,phi] = EnKF(lambda, y1,y2,Nsample,u0,Gamma,Gamma_inv,K,Nmax_it)

%the EnKF procedure for scalar lambda
% y1,y2 K times nObs
% Nsamples=J
% u0 is d times J (here row vector since d=1)
% Gamma is K times K matrix

eps=0.001;     %tolerance
N=0;          
Dt=0.05; 
un=u0;
unp1=un;
unj=u0;
u0_bar = mean(u0, 2);
dy = y(y1,y2,lambda)-G_scalar(lambda,u0_bar,K);
phi0=0.5 * dy' * Gamma_inv * dy; 

phi = eps * phi0 + 1;
while (N<Nmax_it)
    
         N=N+1;
 
         CU=C(un,lambda,Nsample,K);  
         
         yk=y(y1,y2,lambda);  %K times nObs (nObs=1) 
         
         
         for i=1:Nsample  
            % this is for computing u_{n+1},which is a d times J vector
             Gi=G_scalar(lambda,un(:,i),K);  
             unp1(:,i)=un(:,i)+Dt*CU*(Gamma_inv)*(yk-Gi);  
         end
         
         unp1_bar = mean(unp1, 2);
         
         dy = y(y1,y2,lambda)-G_scalar(lambda,unp1_bar,K);
         phi=0.5 * dy' * Gamma_inv * dy;  
         
         un=unp1;
         unj=[unj ; un];  

end

mea=mean(unp1, 2);

end

