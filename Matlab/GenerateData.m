function [Individuals_extra,z_HH_extra_size, HHData_extra_size, number_of_generation] = ...
    GenerateData(hh_size,lambda, w, ...
    phi,pi, d, cum_number_of_generation,possiblehhcount,howmany)
    
    Individuals_extra = [];
    z_HH_extra_size = [];
    HHData_extra_size = [];
    number_of_generation = 0;
    n_possible_household = 0;
    p = length(d);
    while (n_possible_household< possiblehhcount)
        number_of_generation = number_of_generation + 1;
        
        %generate a batch of 10K household
        rn = rand(howmany,hh_size *(1+p) + 2); %the exact number of random numbers needed
        data_to_check = samplehouseholds(phi,w, pi, d, ...
            lambda{1},lambda{2},number_of_generation+cum_number_of_generation, howmany,hh_size,rn);

        %impossible household
        [Households, Index] = checkconstraints(data_to_check);
        possible = howmany - size(Households,1);
        if ((n_possible_household + possible) >= possiblehhcount) 
            Households = Households(1:Index(possiblehhcount - n_possible_household),:);
            
        end
        n_possible_household = n_possible_household + possible;
         
        t1 = Households(:,1:hh_size * 8)';
        t2 = Households(:,hh_size * 8 +1)';
        t3 = Households(:,hh_size * 8 +1 + (1:hh_size))';
        
        
        Individuals = [reshape(t1,8,[])' reshape(repmat(t2,hh_size,1),[],1)...
            reshape(t3,[],1)];
    
        Individuals_extra = [Individuals_extra; Individuals];
        z_HH_extra_size = [z_HH_extra_size  t2];
        HHData_extra_size = [HHData_extra_size Households(:,8)'];
        
    end
    
    number_of_generation = number_of_generation + cum_number_of_generation;
    
    
    
end
