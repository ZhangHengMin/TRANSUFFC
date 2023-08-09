function [acc,NMI,idx] = clustering_liu(Z,gnd)
%
K = max(gnd);
%post processing
[U,S,V] = svd(Z,'econ');
S = diag(S);
r = sum(S>1e-4*S(1));
U = U(:,1:r);S = S(1:r);
U = U*diag(sqrt(S));
U = normr(U);
L = (U*U').^4;

% spectral clustering
D = diag(1./sqrt(sum(L,2)));
L = D*L*D;
[U,S,V] = svd(L);
V = U(:,1:K);
V = D*V;
idx = kmeans(V,K,'emptyaction','singleton','replicates',3,'display','off');
% n = size(V,1);
% M = zeros(K,K,20);
% rand('state',123456789);
% for i=1:size(M,3)
%     inds = false(n,1);
%     while sum(inds)<K
%         j = ceil(rand()*n);
%         inds(j) = true;
%     end
%     M(:,:,i) = V(inds,:);
% end
% idx = kmeans(V,K,'emptyaction','singleton','start',M,'display','off');
% idx 640x1 gnd 1x640
if K>10
    acc  = 100*compacc(idx,gnd);
else
    acc =  100*(1 - missclassGroups(idx,gnd,K)/length(idx));
end
NMI =  100*MutualInfo(gnd,idx);