function stBlockProcess = vocoder_process_pitchmarks(stBlockProcess)
% the vocoder process function by sg
% this function can be used in the ola_framwork

% if the struckt is empty (done by the ola_framework if stBlockProcess was not
% defined in the input args) the construckt the needed struckture:
if isempty(stBlockProcess)
    stBlockProcess = vocoder_init;
end

w = stBlockProcess.w;
f0vec = stBlockProcess.f0vec;
thisInBlock = w.*stBlockProcess.thisInputBlock;
blockLen = size(thisInBlock,1);
[a,k] = lpc_coeff_levinsons_algo(thisInBlock,stBlockProcess.M);
errSig = filter(a,1,thisInBlock);
errSigPower = sum(errSig.^2);
[f0,~,isVoiced] = easy_fundamental_freq_estimate(...
    thisInBlock,stBlockProcess.fs);
f0vec(1) = f0;
f0vec(2:3) = f0vec(1:2);
f0 = median(f0vec);
if  ~isVoiced
    f0 = 0;
end
stBlockProcess.f0_tracker(stBlockProcess.blockNumber)=f0;
if isVoiced
    dkOld = stBlockProcess.dkPulsetrain;
    mPhiOld = stBlockProcess.mPhiPulsetrain;
    stBlockProcess.dkPulsetrain = round(stBlockProcess.fs/f0);
%     stBlockProcess.mPhiPulsetrain = stBlockProcess.dkPulsetrain;    
    stBlockProcess.mPhiPulsetrain = abs(mPhiOld + (dkOld*round(stBlockProcess.nShift/dkOld))-stBlockProcess.nShift);
    if stBlockProcess.mPhiPulsetrain > stBlockProcess.nShift
        stBlockProcess.mPhiPulsetrain = stBlockProcess.mPhiPulsetrain-stBlockProcess.nShift;
    elseif stBlockProcess.mPhiPulsetrain < 1
        stBlockProcess.mPhiPulsetrain = 1;%stBlockProcess.dkPulsetrain;
    end 
%     s = stBlockProcess.s(stBlockProcess.idx);
    s = zeros(blockLen,1);
    s(pitchmarker(thisInBlock,stBlockProcess.fs))=1;
    s = filter([1 .8],1,s); 
else
    s = randn(blockLen,size(thisInBlock,2));
end
% s=errSig; % this would be like the original signal
s = s./max(abs(s));
% thisBlock = filter(2*errSigPower,a,s);
thisBlock = latcfilt(k,errSigPower,s);
thisBlock = w.*thisBlock;
stBlockProcess.thisOutputBlock = thisBlock; 
end