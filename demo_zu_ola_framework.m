%% ideen zu overlap add framework:
% gwünscht ist eine funktion, in die ein procces-funktion function handle
% reingegeben wird.
% außerdem input und output routine soll auch vraiabel sein: zuerst aber
% nur: in variable speichern und input ist auch variable.
clear;
close all;

M = 18;
lenWin = 32e-3;
poverlap = 2/4;
[input,fs] = wavread('Mann.wav');
input = input(1:4*fs,1);
% input = filter([1,-.33],1,input);
input = input(:,1)./max(abs(input(:,1)));
ds = 4;
input = downsample(input,ds);
fs = fs/ds;
nWin = round(lenWin*fs); 
window = hann(nWin+mod(nWin,2),'periodic');
noverlap = round(poverlap*nWin); noverlap = mod(noverlap,2)+noverlap;

% example for process function: direct putting through blocks
stBlockProcess = vocoder_init(M,window,noverlap,fs);
idpm = pitchmarker(input,fs);
s = zeros(size(input,1)+2*length(window)+mod(length(input),length(window)),...
    size(input,2)); 
s(idpm+length(window)+1)=1;
% s(idpm+length(window)+2)=-0.125;
% [b,a]=butter(2,2000/fs,'low');
% s = filter(b,a,s);
s = filter([1 .8],1,s); 

stBlockProcess.s = s;
tic
stOutput = ola_framework(input,@vocoder_process,stBlockProcess); % ola framework 
stOutputPitchmarker = ola_framework(input,@vocoder_process_pitchmarks,stBlockProcess); % ola framework 
toc
figure; plot(stOutput.stBlockProcess.f0_tracker);  
figure; plot(([input,stOutput.output]))

err = sum((input-stOutput.output).^2)
err = sum((input-stOutputPitchmarker.output).^2)
% soundsc(stOutput.output,fs)
soundsc(stOutputPitchmarker.output,fs)