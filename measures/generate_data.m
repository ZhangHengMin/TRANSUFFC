function [X,cids] = generate_data(n,d,D,s,pp,rr)

% n = 20;    % samples in each subspace
% d = 5;     % rank
% D = 200;    % dimension
% s = 15 ;    % number of subspaces


[U,S,V] = svd(rand(D));
cids = [];
U1 = U(:,1:d);
X = U1*rand(d,n);
cids = [cids,ones(1,n)];

for i=2:s
    R = orth(rand(D));
    U1 = R*U1;
    X = [X,U1*rand(d,n)];
    cids = [cids,i*ones(1,n)];
end
nX = size(X,2);
norm_x = sqrt(sum(X.^2,1));
norm_x = repmat(norm_x,D,1);
gn = norm_x.*randn(D,nX);
inds = rand(1,nX)<=pp;
X(:,inds) = X(:,inds) + rr*gn(:,inds);