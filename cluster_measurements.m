function [clusteracc, clusternmi, idx] = cluster_measurements(Z,gnd, ClusterWay)

switch (ClusterWay)
    case 'LiuCai13'
        K = max(gnd);

        for ii = 1 : size(Z,2)
            Z(:,ii) = Z(:,ii) / max(abs(Z(:,ii))) ;
        end
 
        Z = Selection(Z , K ) ;
        Z = ( abs(Z) + abs(Z') ) / 2 ;
 

        [clusteracc, clusternmi, idx] = clustering_liu(Z,gnd);

    case 'Canlucan'
        K = max(gnd);
        for ii = 1 : size(Z,2)
            Z(:,ii) = Z(:,ii) / max(abs(Z(:,ii))) ;
        end

        L = Selection(Z, K) ;
        Z = ( abs(L) + abs(L') ) / 2 ;
        idx = clu_ncut(Z,K) ;
        clusteracc = 100*compacc(idx',gnd);
        clusternmi = 100*compute_nmi(gnd, idx);
end

