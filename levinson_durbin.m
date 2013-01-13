function a = levinson_durbin(r)
% solves the yule walker equation by levinson durbin algorithm
% (that uses properties of toeplitz structured correaltion matrix)
%
% usage:  
%       a = levinson_durbin(r)
%
% input:  
%       r ... (auto)correlation vector (L+1 x 1), where L is the order of
%               the prediction filter.
%
% output: 
%       a ... linear prediction filter coefficients (L+1 x 1)
%
% author: siegfried gündert at uni-oldenburg.de
%
% history: 
%       10.12.2012 <first buildt sg>
%       11.12.2012 <improovement and corrections by SG>
%
% literature: 
%       Handbook of Speech Processing, Chapter 7 "Linear Prediction"

L = length(r)-1;   % Order of the Filter 
a = zeros(L+1,1);  % initialize filter coefficient vector
kappa = r(2);      % initial kappa
MSE   = r(1);      % Mean Squared Error r[0]
K     = -kappa/MSE;    
MSE   = MSE+kappa*K;
a(1)  = K;

for ll = 2:L
    kappa = r(2:ll)'*flipud(a(1:ll-1)) + r(ll+1);
    K = -kappa/MSE;
    a(1:ll) = [a(1:ll-1) ; 0]  +  K*[flipud(a(1:ll-1)) ; 1];    
    MSE = MSE + kappa*K;
end

a(2:L+1) = a(1:L); 
a(1) = 1;