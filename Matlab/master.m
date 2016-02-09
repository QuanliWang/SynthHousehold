init;

%%%%---Blocked sampler---%%
for i = 1:nrun  
    %%
    UpdateParameters;
    cum_number_of_generation = 0;
     %% size2
    hh_size = 2;
    [size2data_to_check_all,size2outcome_to_save, cum_number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, cum_number_of_generation,L,ACS_count, 150000);        
    [Individuals2,HouseHolds2] = GetExtraIndividualAndHousehold(size2data_to_check_all,size2outcome_to_save,hh_size,ACS_count);
    %% size3
    hh_size = 3;
    [size3data_to_check_all,size3outcome_to_save, cum_number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, cum_number_of_generation,L,ACS_count,500000);
    [Individuals3,HouseHolds3] = GetExtraIndividualAndHousehold(size3data_to_check_all,size3outcome_to_save,hh_size,ACS_count);
    %% size4
    hh_size = 4;
    [size4data_to_check_all,size4outcome_to_save, cum_number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, cum_number_of_generation,L,ACS_count,5000000);
    [Individuals4,HouseHolds4] = GetExtraIndividualAndHousehold(size4data_to_check_all,size4outcome_to_save,hh_size,ACS_count);
   
     
    %%    
    %addjust household index, we should remove this code later and
    %automatically track that.
    hh_size_new= zeros(4,1);
    hh_size_new(2) = size(HouseHolds2,1); hh_size_new(3) = size(HouseHolds3,1);hh_size_new(4) = size(HouseHolds4,1);
    
    hh_index = [];cumsize = 0;
    for hh = 2:4
        hh_index = [hh_index; cumsize + reshape(repmat(1:hh_size_new(hh),hh,1),[],1)];
        cumsize = cumsize + hh_size_new(hh);
    end    
    trueextras = [Individuals2; Individuals3; Individuals4];   
    trueextras(:,1) = 10000 + hh_index;
    data_full = [origdata(:,1:8);trueextras(:,1:8)];
    clear hh_index hh cumsize
    trueextras_ownership = [HouseHolds2(:,8); HouseHolds3(:,8); HouseHolds4(:,8)];
    z_HH_new = [HouseHolds2(:,2 * 8 + 1); HouseHolds3(:,3 * 8 + 1); HouseHolds4(:,4*8+1)];
    HHdata1 = [HHdataorig(:,1);trueextras_ownership];
    HHdata2 = [HHdataorig(:,2)-1;ones(hh_size_new(2),1);2*ones(hh_size_new(3),1);3*ones(hh_size_new(4),1)];
    clear trueextras_ownership
    
    
    %%   
    postsave;
    
end

toc    
