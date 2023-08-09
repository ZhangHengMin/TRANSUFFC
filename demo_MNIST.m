clear all;
close all;

addpath('measures');

ii = 1;


dataset='MNISTS_28x28'; %784x2000
data=load([dataset,'.mat']);
X = double(data.fea');
X = X./repmat(sqrt(sum(X.^2)),[size(X,1) 1]);%m*N
gnd = data.gnd';
gnd = double(gnd);


% choices of objective function
ourfuns = {'SpULRR1', 'SpULRR12', 'SpULRR23'};
for ifun =  1 : length(ourfuns)
    funname = ourfuns{ifun};
    disp([' funname= ' num2str(funname)]);
    % lambda
    La = [1e-4, 1e-3, 1e-2, 1e-1 1.0];
    for k = 1  :length(La)
        lambda = La(k);
        %disp([' lambda = ' num2str(lambda)]);
        % rank
        rankParas = [20 30 40 50 60];
        for inparas = 1   : length(rankParas)
            ranknum = rankParas(inparas);% lambda


            %% with optimal mean
            switch(funname)
                case 'SpULRR1'
                    tic;
                    [Z, iter, stopC1, Time1] = SpULRR1(X, lambda, ranknum, 1e-5);
                    timecost = toc;
                case 'SpULRR12'
                    tic;
                    [Z, iter, stopC2, Time2] = SpULRR12(X, lambda, ranknum, 1e-5);
                    timecost = toc;
                case 'SpULRR23'
                    tic;
                    [Z, iter, stopC3, Time3] = SpULRR23(X, lambda, ranknum, 1e-5);
                    timecost = toc;


            end

            %% clustering results
            Clusterfuns = {'LiuCai13', 'Canlucan'};
            ClusterWay = Clusterfuns{ii};
            [clusteracc, clusternmi] = cluster_measurements(Z,gnd, ClusterWay);

            disp([' ClusterWay = ' num2str(ClusterWay), ' numrank = ' num2str(ranknum), ' method = ' num2str(funname),...
                '  acc= ' num2str(clusteracc), '  nmi= ' num2str(clusternmi), '  time= ' num2str(timecost), '  iter= ' num2str(iter)]);


        end
 
    end
end



