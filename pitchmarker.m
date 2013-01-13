function [m]=pitchmarker(in,Fs)

% Usage: PSOLA(In,m,Fs,alpha,beta)
% PSOLA:    Time Stretches Audio Input by Faktor X
% -------------------------------------------------------------------------
% 
% 
%   input:   ---------
%            
%            
%  output:   ---------
%            
%
% Author:  Tobias Bruns Mat.Nr.: 5019011  
% Version: April 17 2010

m = zeros(length(in),1);% Allocation of Pitchmark Vevtor
f0=150; % Hz
w0 = 2*pi*f0/Fs; %1/s
Q = 10; % Guete
alpha_param = sin(w0)/(2*Q); %Filter Parameter alpha

% Calculate Biquad Coefficiants

% Bandpass
b0 =   Q*alpha_param;
b1 =   0;
b2 =   -Q*alpha_param;
a0 =   1 + alpha_param;
a1 =  -2*cos(w0);
a2 =   1 - alpha_param;

% Tiefpass
% b0 =  (1 - cos(w0))/2;
% b1 =   1 - cos(w0);
% b2 =  (1 - cos(w0))/2;
% a0 =   1 + alpha_param;
% a1 =  -2*cos(w0);
% a2 =   1 - alpha_param;

% Filter of Processed Audio Block
ProcessBlock = filter([b0 b1 b2],[a0 a1 a2],in);
% soundsc(ProcessBlock,Fs);
% figure(1);pwelch(ProcessBlock)
% figure(1);freqz([b0 b1 b2], [a0 a1 a2])

% Sign to find Zero Crossings
ProcessBlock = sign(ProcessBlock);
ProcessBlock(ProcessBlock==0) = 1;
m = ProcessBlock(1:end-1) - ProcessBlock(2:end);
m(m>0) = 0;
m = find(abs(m) > 0);

% experiment: filter too short pitchmarks, maybe also filter too long marks
pit = diff(m);
idx = find(pit>median(pit)/2);
m   = m(idx+1);

% Generate Pitchmarks
% for kk = 2: length(ProcessBlock)
%     if ProcessBlock(kk)>=ProcessBlock(kk-1)
%         m(kk)=1;
%     end
% end

% figure(2);
% hold on;
% plot(m);
% plot(in,'r');
% hold off;
