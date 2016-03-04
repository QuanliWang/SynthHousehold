function [z_Households_individial_extra,HouseHolds_extra,phicount,number_of_generation] = GenerateData2(hh_size,lambda, w, ...
    phi,pi, d, p, cum_number_of_generation,L,ACS_count)
    maxd = size(phi,1)/ p;
    K= size(w,1);
    phicount=zeros(K*L,maxd,p);      % phi counts
    
    z_Households_individial_extra = [];
    HouseHolds_extra = [];
    number_of_generation = 0;
    n_possible_household = 0;
    while (n_possible_household< ACS_count(hh_size-1))
        number_of_generation = number_of_generation + 1;
        
        %generate a batch of 10K household
        data_to_check = GenerateData2Check(hh_size,lambda, w, ...
            phi,pi, d, p, number_of_generation+cum_number_of_generation, L);
        
        % transpose the data to check
        data_to_checkT = data_to_check(:,1:(hh_size*8))';   
        outcome = checkingconstraints(data_to_checkT,10000,hh_size);
        possible = sum(outcome);
        
        %impossible household
        HouseHolds = data_to_check(outcome==0,:);

        if ((n_possible_household + possible) >= ACS_count(hh_size-1)) 
            %possible household
            possiblehh = data_to_check(outcome==1,1);
            %find the number ACS_count(hh_size-1)'th possible household 
            exceeding_index = possiblehh(ACS_count(hh_size-1) - n_possible_household,1);
            %Get all the impossible households we need
            HouseHolds = HouseHolds(HouseHolds(:,1) < exceeding_index,:);
        end
        n_possible_household = n_possible_household + possible;
         
        Individuals = reshape(HouseHolds(:,1:hh_size * 8)',8,[]);
        z_hh_individual1 = reshape(repmat(HouseHolds(:,hh_size * 8 +1)',hh_size,1),[],1);
        z_hh_individual2 = reshape(HouseHolds(:,hh_size * 8 +1 + (1:hh_size))',[],1);
       
        % add code to do partial counts
        groupIndex = L*(z_hh_individual1-1)+z_hh_individual2;
        dataShift = 2; 
        
        for j = 1:p
            phicount(:,1:d(j),j) = phicount(:,1:d(j),j) + ...
                groupcount(groupIndex,Individuals(dataShift+j,:),K*L, d(j));
        end
    
        HouseHolds_extra = [HouseHolds_extra;HouseHolds];
        z_Households_individial_extra = [z_Households_individial_extra;...
                    [z_hh_individual1 z_hh_individual2]];
        
    end
    
    number_of_generation = number_of_generation + cum_number_of_generation;
end
