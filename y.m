function [meas] = y(y1,y2,lambda)

% Compute the function y=lambda*y1+(1-lambda)*y2 
% y_i is in k \times nObs
% k is the dimension of the observation space 
% nObs number of observations 

if isequal(size(y1),size(y2))

    meas=lambda*y1+(1-lambda)*y2;
else
    meas = NaN;
    MException('MATLAB:EnKF:y:InputError', ' Input dimensions for multiobjective y do not agree ')

end

end

