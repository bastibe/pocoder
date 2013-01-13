function [f0 mm,isVoiced] = easy_fundamental_freq_estimate(thisInBlock,fs)

flag_verbose = 0;
g_noisyness = 1.9;

f0_erw = [80,230];
f0lag  = round(fs./f0_erw);
rxx = xcorr( thisInBlock-mean(thisInBlock) , f0lag(1) );
[mm,ii] = max( rxx((f0lag(2):f0lag(1))+round(length(rxx)/2)) );
iiMax   = f0lag(2)  + ii - 1;
T0 = iiMax/fs;
f0 = 1/T0;
isVoiced = sum(diff(sign(diff(rxx)))<0)<(abs(diff(f0lag)/g_noisyness));  

if flag_verbose
    figure(1); 
    plot( rxx );
    hold on; plot( iiMax+round(length(rxx)/2 ) , rxx( iiMax+round(length(rxx)/2) ),'ro')
    hold off;
    disp(sprintf('Fundamental frequency current block: %.2f',f0))
end