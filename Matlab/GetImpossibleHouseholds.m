function [HHdata_individual_extra,z_HH_extra,IndividualData_extra,HHdata_extra,hh_size_new] = ...
    GetImpossibleHouseholds(d,ACS_count,lambda,w,phi,pi)
    L = size(w,2); %get 
    p = length(d);
    %%
    Individuals = cell(4,1);
    HouseHolds =cell(4,1);
    cum_number_of_generation = 0;
    for hh_size = 2:4
        [Individuals{hh_size},HouseHolds{hh_size}, cum_number_of_generation] = ...
            GenerateData2Check_all(hh_size,lambda, w, ...
            phi,pi, d, p, cum_number_of_generation,L,ACS_count); 
    end
    
    %%    
    cumsize = 0;
    hh_size_new= zeros(3,1);
    hh_index = [];
    ImpossibleIndividuals = [];
    z_HH_extra = [];
    HHdata_extra = [];
    for hh = 2:4
        hh_size_new(hh-1) = size(HouseHolds{hh},1); 
        hh_index = [hh_index; cumsize + reshape(repmat(1:hh_size_new(hh-1),hh,1),[],1)];
        current_household = HouseHolds{hh};
        cumsize = cumsize + hh_size_new(hh-1);
        ImpossibleIndividuals = [ImpossibleIndividuals;Individuals{hh}];
        z_HH_extra = [z_HH_extra; current_household(:,hh * 8 + 1)];
        
        HHdata_extra = [HHdata_extra ; [current_household(:,8) ones(hh_size_new(hh-1),1) * (hh -1)]];
    end 
    ImpossibleIndividuals(:,1) = 10000 + hh_index;
    IndividualData_extra = ImpossibleIndividuals(:,1:8);
    HHdata_individual_extra = ImpossibleIndividuals(:,9:10);
end

