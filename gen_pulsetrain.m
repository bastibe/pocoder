function d = gen_pulsetrain(N,mPhi,dk)
d = zeros(N,1);
% M = round((N/dk)/2);
d(mPhi:dk:end) = 1; 
d = filter([1 .8],1,d); 
