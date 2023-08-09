function [Z,iter,  stopC,Time ]=SpULRR1(X, lambda, ranks, mu, tol, maxIter,term_obj)
% 2022.4.6 update
% min lambda*1/2(\|U\|_F^2+\|V\|_F^2)+|E|_2,p
% s.t. Z=UV^T,X = XZ+b1^T+E
% inputs:
%        X -- m*n data matrix, m is the data dimension, and n is the number
%             of data vectors.
%% Default parameters
[ncr(1,1),ncr(1,2)] = size(X);
max_beta = 1e10;
p = 1;
rho = 1.5;
if nargin < 7
    term_obj = 1;
end
if nargin < 6
    maxIter = 500;
end
if nargin < 5
    tol = 1e-8;
end
if nargin < 4
    mu = 1e-5; % turnable 1e-6,..,1e0
end
if nargin < 3
    ranks = round(1.2*rank_estimation(D));
end
if nargin < 2
    lambda = sqrt(max(ncr));
end
%% Initializing optimization variables
% U=cell(1,2);
% for ii = 1:2
%     U{ii} = rand(ncr(1, 2), ranks);
% end
for ii = 1:2
    if ii == 1
        M{ii} = rand(size(X, ii), ranks);              % fix
        %M{ii} = randn(size(D, ii), ranks);
        [U{ii}, aa1] = qr(X'*M{ii}, 0);
    else
        M{ii} = rand(size(X, 3-ii), ranks);            % fix
        %M{ii} = randn(size(D, 3-ii), ranks);
        [U{ii}, aa2] = qr(X'*M{ii}, 0);
    end
end
Z  = U{1}*U{2}'; 
E = zeros(ncr(1,1), ncr(1,2));
b = zeros(ncr(1,1),1);
W3 = eye(ncr(1,2));

%% Advance compute
XXEye = X'*X + eye(ncr(1,2));
oneV = inv(ones(1,ncr(1,2))*ones(ncr(1,2),1));
%% Start main loop
flagTime=tic;
for iter=1:maxIter
    beta=1/mu;Lbeta=lambda*beta;
    %-----------upadat U{1}=U, U{2}=V----------------%
    u22=U{2}'*U{2};
    U{1}=Z*U{2}/(Lbeta*eye(ranks)+u22);%
    U{2} = Z'*U{1}/(Lbeta*eye(ranks)+U{1}'*U{1});
    %-------------update Z----------------%
    M1=U{1}*U{2}';
    XE=X-E;
    M2=XE-b*ones(1,ncr(1,2));
    Z=XXEye\(M1+X'*M2);
    %---------------update b----------------%
    XZ=X*Z;
    M3=XE-XZ;
    b=M3*ones(ncr(1,2),1)*oneV;
    %---------------update E----------------%
    M4=X-XZ-b*ones(1,ncr(1,2));
    E=M4/(2*beta*W3+eye(ncr(1,2)));
    %---------------update W3---------------%
    Ei1 = sqrt(sum(E.*E,1)+eps)';
    w = (p/2)./(Ei1.^(2-p)); %d1 = 0.5*p./(Xi1.^(p/2));
    W3 = spdiags(w,0,ncr(1,2),ncr(1,2));
    
    %---------------compute converge---------------%
    leq1=M4-E;
    stopC(iter) = max(max(abs(leq1)));

    Time(iter)=toc(flagTime);
    %----------if converge-------%
    if stopC(iter)<tol
        break
    else
        mu = min(mu * rho, max_beta);
    end
end