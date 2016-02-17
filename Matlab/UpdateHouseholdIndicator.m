function [z_HH_Individuals,z_HH] = ...
    UpdateHouseholdIndicator(n,n_individuals,maxd,p,dataT,SS,...
    K,L,phi,w,pi,HHdata,lambda)
    z_HH= zeros(n,1);                  %cluster indicators for household; K
    z_HH_prob = samplezHHwithHHnewv1_2HHvar(phi,dataT,w,pi,SS',...
      K,L,p,maxd,n,HHdata(:,1),lambda{1},HHdata(:,2),lambda{2});
    for h = 1:n
       zupdateprob_h = z_HH_prob(K*(h-1)+1:K*h);
       z_HH(h) = randomsample(zupdateprob_h,rand);
    end
    
    z_HH_Individuals = zeros(n_individuals,1);
    cumsumSS = [0 cumsum(SS)];
    for h = 1:n
       z_HH_Individuals(cumsumSS(h)+1:cumsumSS(h+1)) = z_HH(h);
    end
  
end

