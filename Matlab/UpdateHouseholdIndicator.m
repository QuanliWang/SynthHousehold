function [z_HH_all, z_HH_Individuals_all] = UpdateHouseholdIndicator(n,n_individuals,maxd,p, K,L,...
    phi,dataT,SS,w,pi,HHdata_all,lambda,ImpossibleIndividuals,z_HH_extra)
    z_HH= zeros(n,1);                  %cluster indicators for household; K
    z_HH_prob = samplezHHwithHHnewv1_2HHvar(phi,dataT,w,pi,SS',K,L,p,maxd,n,HHdata_all(1:n,1),lambda{1},HHdata_all(1:n,2),lambda{2});
    for h = 1:n
       zupdateprob_h = z_HH_prob(K*(h-1)+1:K*h);
       z_HH(h) = randomsample(zupdateprob_h,rand);
    end
    z_HH_all = [z_HH;z_HH_extra];
    
    z_HH_Individuals = zeros(n_individuals,1);
    cumsumSS = [0 cumsum(SS)];
    for h = 1:n
       z_HH_Individuals(cumsumSS(h)+1:cumsumSS(h+1)) = z_HH(h);
    end
   
    if isempty(ImpossibleIndividuals)
        z_HH_Individuals_all = z_HH_Individuals;
    else
        z_HH_Individuals_all = [z_HH_Individuals;ImpossibleIndividuals(:,9)];
    end
end

