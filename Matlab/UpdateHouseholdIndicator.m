function [z_HH_all, z_HH_Individuals_all,newphi] = UpdateHouseholdIndicator(n,n_individuals,maxd,p, K,L,...
    phi,dataT,SS,w,pi,HHdata1_all,lambda1,HHdata2_all,lambda2,ImpossibleIndividuals,z_HH_extra)
    z_HH= zeros(n,1);                  %cluster indicators for household; K
   
    %%% to work in C++, consider make phi into a 2-d matrix maxd*p-by-K*L
    newphi = zeros(maxd*p,K*L);
    for e = 1:K*L
       newphimatrix = reshape(phi(e,:,:),p,maxd);
       for j = 1:p
           newphi(maxd*(j-1)+1:maxd*j,e) = newphimatrix(j,:);
       end
    end
    z_HH_prob = samplezHHwithHHnewv1_2HHvar(newphi,dataT,w,pi,SS',K,L,p,maxd,n,HHdata1_all(1:n),lambda1,HHdata2_all(1:n),lambda2);
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

    
    disp('zHH updated');
end
