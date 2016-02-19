function [Individuals_all,HouseHolds_all, number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda, w, ...
    phi,pi, d, p, cum_number_of_generation,L,ACS_count,newxtras)
    
    Individuals_all = [];
    HouseHolds_all = [];
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
         
        t1 = HouseHolds(:,1:hh_size * 8)';
        t2 = HouseHolds(:,hh_size * 8 +1)';
        t3 = HouseHolds(:,hh_size * 8 +1 + (1:hh_size))';
        
        
        Individuals = [reshape(t1,8,[])' reshape(repmat(t2,hh_size,1),[],1)...
            reshape(t3,[],1)];
    
        Individuals_all = [Individuals_all; Individuals];
        HouseHolds_all = [HouseHolds_all;HouseHolds];
        
    end
    
    number_of_generation = number_of_generation + cum_number_of_generation;
    
    
    
end
