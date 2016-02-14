function [hh_size_new,ImpossibleIndividuals,data_full_all,z_HH_extra,HHdata1_all,HHdata2_all] = ...
    GetImpossibleHouseholds(lambda1,lambda2,w,phi,d,p,L, ACS_count,pi,...
    origdata,HHdataorig)
    
    %%
    list_size = [0 150000 500000 5000000];
    Individuals = cell(4,1);
    HouseHolds =cell(4,1);
    cum_number_of_generation = 0;
    for hh_size = 2:4
        [data_to_check_all,outcome_to_save, cum_number_of_generation] = ...
            GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
            phi,pi, d, p, cum_number_of_generation,L,ACS_count, list_size(hh_size));        
        [Individuals{hh_size},HouseHolds{hh_size}] = GetExtraIndividualAndHousehold(data_to_check_all,...
            outcome_to_save,hh_size,ACS_count);
    end
    
    %%    
    cumsize = 0;
    hh_size_new= zeros(4,1);
    hh_index = [];
    ImpossibleIndividuals = [];
    z_HH_extra = [];
    HHdata1_all = HHdataorig(:,1);
    HHdata2_all = HHdataorig(:,2)-1;
    for hh = 2:4
        hh_size_new(hh) = size(HouseHolds{hh},1); 
        hh_index = [hh_index; cumsize + reshape(repmat(1:hh_size_new(hh),hh,1),[],1)];
        current_household = HouseHolds{hh};
        cumsize = cumsize + hh_size_new(hh);
        ImpossibleIndividuals = [ImpossibleIndividuals;Individuals{hh}];
        z_HH_extra = [z_HH_extra; current_household(:,hh * 8 + 1)];
        HHdata1_all = [HHdata1_all; current_household(:,8)];
        HHdata2_all = [HHdata2_all; ones(hh_size_new(hh),1) * (hh -1)];
    end 
    ImpossibleIndividuals(:,1) = 10000 + hh_index;
    data_full_all = [origdata(:,1:8);ImpossibleIndividuals(:,1:8)];
    disp('synthesis updated');
end
