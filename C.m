function [c] = C(u,lambda,Nsample,K)

% Compute the covariance matrix

% u is a matrix d \times J
% d=1 dimension of the settings of the machine
% J=Nsample is the number of samples
% lambda is a scalar
% K number of the observations (dimension of y)


ubar=mean(u,2);  
Gbar=zeros(K,1);
for i=1:Nsample
 Gbar=Gbar+G_scalar(lambda,u(:,i),K);  
end
Gbar=Gbar/Nsample;
c=zeros(size(u,1),K);   %d \times K

for i=1:Nsample
c=c+(u(:,i)-ubar)*(G_scalar(lambda,u(:,i),K)-Gbar)';

end

c=c/Nsample;

end
