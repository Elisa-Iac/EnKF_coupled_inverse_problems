function g2 = G2_scalar(u, K)
   % maps u in R^1 to y in R^K
   %Examples
   g2 = (u+0.5)^2 * ones(K, 1); %PB1
   %g2=1 - exp(-(u+1).^2);  %PB3
end