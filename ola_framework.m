function stOutput = ola_framework(input,fun_process,stProcess)
% overlap add framework that gets a block process function as input
%
% usage: output = ola_framework(input,fun_process,window,overlap)
%
% input: 
%              input ... is the whole input signal to be processed (N x nChan)
%        fun_process ... function handle to the process function:
%                        example:  fun_process = @(x_in) x_in;
%                                  would reac all samples through tothe
%                                  output.
%             window ... can either be an integer denoting the length of an
%                        preiodic hann window or a vector contaning the 
%                        window samples.
%          stProcess ... Structure for the process function. Normally you
%          dont need to define it here if you initialize your struckt in
%          your process function in the case it is empty. but if you use
%          the ola_framework in another blocking loop you mabe need to give
%          through the states from one block to another.
% 
% output:
%             output ... is the whole output signal (N x nChan)
%
% author: siegfried guendert {at} uni-oldenburg.de
%
% history: 20.12.2012 <first buildt, sg>
%% WINDOW NOCH RAUSNEHMEN UND IN PROCESS VERSCHIEBEN!!!!
%
% examples:
%
%
% parse input
if nargin<2, error('nargin:stft: check input arguments'); end

% check window: if it is a scalar value, the generate a periodic "von hann"
% window with the length [window]:
% if numel(window) == 1
%     nW = round(window);
%     window = hann(nW+mod(nW,2),'periodic');
% end

blockLen = stProcess.blockLen;
if mod(blockLen,2)
    warning('ola_framework:windowLength',...
        ['blockLen: the length of your window should be an even number ',...
        'to prevent amplidude modulation by the window. By the way ',...
        'your window should be ''periodic''. Alternatively use a ',...
        'scalar value as window to define the window length, the rest,'...
        'is done automatically']);
end
% if nargin == 3
%     nOverlap = round(blockLen/2);
% elseif overlap<1
%     nOverlap = round(overlap*blockLen);
% else
%     nOverlap = overlap;        
% end
nOverlap = stProcess.nOverlap;
if nargin <3
    stProcess = [];
end
%     
% number of samples the sliding window is shifted
nShift = round(blockLen-nOverlap);
stProcess.nShift = nShift;
nSamples = size(input,1);
nChan = size(input,2);

% zeropad the signal so we can reconstruct the whole signal in the end:
x = [ zeros(blockLen,nChan);...
      input;...
      zeros(blockLen+mod(nSamples,blockLen),nChan)];
y = zeros(size(x));  

% repeat windows nChan times:
% window = repmat(window(:),1,nChan);  

% calculate number of blocks to be computed:
sizX = size(x);
nBlocks = floor((sizX(1)-nOverlap)/nShift);

% do block process:
for kk = 1:nBlocks        
    %% reading stage:    
    idx = (1:blockLen) + (kk-1)*nShift;   % current block index vector   
    stProcess.idx = idx;    
    thisInputBlock = x(idx,:);            % current  block:
    stProcess.thisInputBlock = thisInputBlock;
    stProcess.blockNumber = kk;    
    %% process stage:
    stProcess = fun_process(stProcess);
    thisOutputBlock = stProcess.thisOutputBlock;
%     figure(101); plot(idx,thisOutputBlock,'color',rand(1,3)); hold on;
    %% adding overlapped blocks stage:
    y(idx,:) = y(idx,:) + thisOutputBlock;    

end

% throw away the first and last zeropadded block:
output = y((blockLen+1):(blockLen+nSamples),:); 

% calibrate the amplitude: % A = sqrt(sum(x.^2)/sum(y.^2)); % empiric
A = nShift/blockLen; % this is the analytic solution. so 
                       % it doesent change influences of the process- 
                       % functions to the amplitude
                     
stOutput.output = output.*A;
stOutput.stBlockProcess=stProcess;