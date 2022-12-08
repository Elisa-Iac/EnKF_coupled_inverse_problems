function [g] = G_scalar(lambda1,u,K)

% Compute the function G(u)=lambda*G1+(1-lambda)*G2
% u is a vector in R^1
% lambda is a scalar
% a vectorized shape is required for the rest of the implementation

g=lambda1*G1_scalar(u, K) + (1-lambda1)*G2_scalar(u, K);

end
