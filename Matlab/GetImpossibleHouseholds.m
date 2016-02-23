function [HHdata_individual_extra,z_HH_extra,IndividualData_extra,HHdata_extra,hh_size_new] = ...
    GetImpossibleHouseholds(d,ACS_count,lambda,w,phi,pi,howmany,n)
    L = size(w,2);
    p = length(d);
    
    %%
    cumsize = 0;
    hh_size_new= zeros(3,1);
    hh_index = [];
    ImpossibleIndividuals = [];
    z_HH_extra = [];
    HHdata_extra = [];
    %%
    cum_number_of_generation = 0;
    for hh_size = 2:4
        [Individuals,z_HH_extra_size, HHData_extra_size, cum_number_of_generation] = ...
            GenerateData(hh_size,lambda, w, ...
            phi,pi, d, cum_number_of_generation,ACS_count(hh_size -1),howmany); 
        
        hh_size_new(hh_size-1) = length(z_HH_extra_size); 
        hh_index = [hh_index; cumsize + reshape(repmat(1:hh_size_new(hh_size-1),hh_size,1),[],1)];
        cumsize = cumsize + hh_size_new(hh_size-1);
        ImpossibleIndividuals = [ImpossibleIndividuals Individuals];
        z_HH_extra = [z_HH_extra z_HH_extra_size];
        HHdata_extra = [HHdata_extra [HHData_extra_size; ones(1,hh_size_new(hh_size-1)) * (hh_size -1)]];
        
    end
    
    %%    
    ImpossibleIndividuals(1,:) = n + hh_index;
    IndividualData_extra = ImpossibleIndividuals(1:8,:);
    HHdata_individual_extra = ImpossibleIndividuals(9:10,:);
end

