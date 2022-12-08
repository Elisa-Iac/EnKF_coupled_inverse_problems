function g1 = G1_scalar(u, K)
    % maps u in R^1 to y in R^K
    % Examples
   
    g1 = (u-0.5)^2 * ones(K, 1); %PB1 %[-1;1] Convex case
    %g1= 1 - exp(-(u-1).^2);  %PB3  %[-4,4]  Concave case
end