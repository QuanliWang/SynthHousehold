function [z_Individual_all] = UpdateIndividualIndicator(n,n_individuals,...
    phi,dataT,w,K,L,p,maxd,z_HH_all,HHserial,z_Individual_extra)
    %% update zmember
    %%% to work in C++
    z_member = zeros(n_individuals,1);             %cluster indicators for household member; L
    z_member_prob = samplezmemberv1(phi,dataT,w,K,L,p,maxd,n_individuals,z_HH_all(1:n),HHserial);

    for m = 1:n_individuals
       zupdateprob_m = z_member_prob(L*(m-1)+1:L*m);
       z_member(m) = randomsample(zupdateprob_m,rand);
    end
    z_Individual_all = [z_member;z_Individual_extra(:,2)];
end

