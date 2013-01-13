function [a,k] = lpc_coeff_levinsons_algo(s,n)
% compute LP coefficients a of s by the levonson durbin algorithm
% 
% usage: a = lpc_coeff(s,n)
% iput: ------------------
%       s... input signal
%       n... Order of the AR Filter
% output:------------------
%       a... lpc filter coefficients
%
% yule walker equations by the levinson durbin algorithm
%
% author: siegfried gündert
% history: first version 04.12.2012
if numel(n)>1,error('lpc_coeff:numel(n): n has to be a scalar value'); end

% cross correlation function with order lag n
r_ss = xcorr(s(:),n);

% with levinson durbin algorithm: use of r[0]...r[n]
[a,k] = levinson_algorithmus(r_ss(n+1:end));