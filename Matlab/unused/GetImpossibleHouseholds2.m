function [z_HHdata_individual_extra,z_HH_extra,HHdata_extra,picount_extra,hh_size_new] = ...
    GetImpossibleHouseholds2(d,ACS_count,lambda,w,phi,pi)
    
    L = size(w,2);
    p = length(d);
    
    %%
    hh_size_new= zeros(3,1);
    z_HHdata_individual_extra = [];
    z_HH_extra = [];
    HHdata_extra = [];
    picount_extra = 0;
    %%
    cum_number_of_generation = 0;
    for hh_size = 2:4
        [z_Individuals,HouseHolds, picount_partial,cum_number_of_generation] = ...
            GenerateData2(hh_size,lambda, w, ...
            phi,pi, d, p, cum_number_of_generation,L,ACS_count); 
        
        hh_size_new(hh_size-1) = size(HouseHolds,1); 
        z_HHdata_individual_extra = [z_HHdata_individual_extra; z_Individuals];
        z_HH_extra = [z_HH_extra; HouseHolds(:,hh_size * 8 + 1)];
        HHdata_extra = [HHdata_extra ; [HouseHolds(:,8) ones(hh_size_new(hh_size-1),1) * (hh_size -1)]];
        picount_extra = picount_partial + picount_extra;
    end
   
end

