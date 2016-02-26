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
        [outcome, Households, Index, synHouseholds, possible] = checkconstraints(data_to_check,possiblehhcount-n_possible_household);
        n_possible_household = n_possible_household + possible;
        
        Individuals = households2individuals(Households);
        
        %the c code is actrually slower
        %Individuals_extra = mergeindividuals(Individuals_extra,Individuals);
        Individuals_extra = [Individuals_extra Individuals];
        
        z_HH_extra_size = [z_HH_extra_size  Households(hh_size * 8 +1,:)];
        HHData_extra_size = [HHData_extra_size Households(8,:)];
        
    end
    
    number_of_generation = number_of_generation + cum_number_of_generation;
    
end
