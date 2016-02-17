function z_Individuals = UpdateIndividualIndicator(p,maxd,dataT,HHserial,...
    K,L,z_HH,phi,w)
    n_individuals = size(dataT,2);
    %% update zmember
    %%% to work in C++
    z_Individuals = zeros(n_individuals,1);             %cluster indicators for household member; L
    z_Individuals_prob = samplezmemberv1(phi,dataT,w,K,L,p,maxd,n_individuals,z_HH,HHserial);

    for m = 1:n_individuals
       zupdateprob_m = z_Individuals_prob(L*(m-1)+1:L*m);
       z_Individuals(m) = randomsample(zupdateprob_m,rand);
    end
end

