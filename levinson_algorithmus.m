function [a,k]=levinson_algorithmus(r)

a = zeros(length(r),1);
k = zeros(length(r),1);
for m = 1:length(r)-1 
    alpha = -((r(m:-1:1).' * [a(m:-1:2);1]));
    mu = -(r(m:-1:1).'*[a(2:m);0]) - r(m+1);
    k(m) = -mu/alpha;     
    a(2:m+1) = [a(2:m);0] + k(m)*[a(m:-1:2);1];
end
a(1) = 1;
