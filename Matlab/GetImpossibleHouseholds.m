function [hh_size_new,trueextras,data_full_all,z_HH_extra,HHdata1_all,HHdata2_all] = ...
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
    
%     %% size2
%     hh_size = 2;
%     [size2data_to_check_all,size2outcome_to_save, cum_number_of_generation] = ...
%     GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
%     phi,pi, d, p, cum_number_of_generation,L,ACS_count, 150000);        
%     [Individuals2,HouseHolds2] = GetExtraIndividualAndHousehold(size2data_to_check_all,size2outcome_to_save,hh_size,ACS_count);
%     %% size3
%     hh_size = 3;
%     [size3data_to_check_all,size3outcome_to_save, cum_number_of_generation] = ...
%     GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
%     phi,pi, d, p, cum_number_of_generation,L,ACS_count,500000);
%     [Individuals3,HouseHolds3] = GetExtraIndividualAndHousehold(size3data_to_check_all,size3outcome_to_save,hh_size,ACS_count);
%     %% size4
%     hh_size = 4;
%     [size4data_to_check_all,size4outcome_to_save, cum_number_of_generation] = ...
%     GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
%     phi,pi, d, p, cum_number_of_generation,L,ACS_count,5000000);
%     [Individuals4,HouseHolds4] = GetExtraIndividualAndHousehold(size4data_to_check_all,size4outcome_to_save,hh_size,ACS_count);
%    
     
    %%    
    cumsize = 0;
    hh_size_new= zeros(4,1);
    hh_index = [];
    trueextras = [];
    z_HH_extra = [];
    HHdata1_all = HHdataorig(:,1);
    HHdata2_all = HHdataorig(:,2)-1;
    for hh = 2:4
        hh_size_new(hh) = size(HouseHolds{hh},1); 
        hh_index = [hh_index; cumsize + reshape(repmat(1:hh_size_new(hh),hh,1),[],1)];
        current_household = HouseHolds{hh};
        cumsize = cumsize + hh_size_new(hh);
        trueextras = [trueextras;Individuals{hh}];
        z_HH_extra = [z_HH_extra; current_household(:,hh * 8 + 1)];
        HHdata1_all = [HHdata1_all; current_household(:,8)];
        HHdata2_all = [HHdata2_all; ones(hh_size_new(hh),1) * (hh -1)];
    end 
    trueextras(:,1) = 10000 + hh_index;
    data_full_all = [origdata(:,1:8);trueextras(:,1:8)];
end

