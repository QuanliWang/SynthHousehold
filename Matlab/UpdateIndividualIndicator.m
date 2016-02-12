function [z_Individual_all] = UpdateIndividualIndicator(n,n_individuals,...
    newphi,dataT,w,K,L,p,maxd,z_HH_all,HHserial,ImpossibleIndividuals)
    %% update zmember
    %%% to work in C++
    z_member = zeros(n_individuals,1);             %cluster indicators for household member; L
    z_member_prob = samplezmemberv1(newphi,dataT,w,K,L,p,maxd,n_individuals,z_HH_all(1:n),HHserial);

    for m = 1:n_individuals
       zupdateprob_m = z_member_prob(L*(m-1)+1:L*m);
       z_member(m) = randomsample(zupdateprob_m,rand);
    end

    if isempty(ImpossibleIndividuals)
       z_Individual_all = z_member;
    else
       z_Individual_all = [z_member;ImpossibleIndividuals(:,10)];
    end

    %toc
    disp('zmember updated');
end

