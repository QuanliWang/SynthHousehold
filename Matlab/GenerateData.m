function [Individuals_extra,HouseHolds_extra, number_of_generation] = ...
    GenerateData(hh_size,lambda, w, ...
    phi,pi, d, cum_number_of_generation,ACS_count,howmany)
    
    Individuals_extra = [];
    HouseHolds_extra = [];
    number_of_generation = 0;
    n_possible_household = 0;
    p = length(d);
    while (n_possible_household< ACS_count(hh_size-1))
        number_of_generation = number_of_generation + 1;
        
        %generate a batch of 10K household
        rn = rand(howmany,hh_size *(1+p) + 2); %the exact number of random numbers needed
        data_to_check = samplehouseholds(phi,w, pi, d, ...
            lambda{1},lambda{2},number_of_generation+cum_number_of_generation, howmany,hh_size,rn);

        outcome = checkingconstraints(data_to_check);
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
    
        Individuals_extra = [Individuals_extra; Individuals];
        HouseHolds_extra = [HouseHolds_extra;HouseHolds];
        
    end
    
    number_of_generation = number_of_generation + cum_number_of_generation;
    
    
    
end
