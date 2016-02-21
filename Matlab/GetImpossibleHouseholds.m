function [HHdata_individual_extra,z_HH_extra,IndividualData_extra,HHdata_extra,hh_size_new] = ...
    GetImpossibleHouseholds(d,ACS_count,lambda,w,phi,pi)
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
        [Individuals,HouseHolds, cum_number_of_generation] = ...
            GenerateData(hh_size,lambda, w, ...
            phi,pi, d, p, cum_number_of_generation,L,ACS_count); 
        
        hh_size_new(hh_size-1) = size(HouseHolds,1); 
        hh_index = [hh_index; cumsize + reshape(repmat(1:hh_size_new(hh_size-1),hh_size,1),[],1)];
        cumsize = cumsize + hh_size_new(hh_size-1);
        ImpossibleIndividuals = [ImpossibleIndividuals;Individuals];
        z_HH_extra = [z_HH_extra; HouseHolds(:,hh_size * 8 + 1)];
        HHdata_extra = [HHdata_extra ; [HouseHolds(:,8) ones(hh_size_new(hh_size-1),1) * (hh_size -1)]];
        
    end
    
    %%    
    ImpossibleIndividuals(:,1) = 10000 + hh_index;
    IndividualData_extra = ImpossibleIndividuals(:,1:8);
    HHdata_individual_extra = ImpossibleIndividuals(:,9:10);
end

