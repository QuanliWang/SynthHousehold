function dummy = master()
init;

%profile on
%%%%---Blocked sampler---%%
for i = 1:nrun    
    UpdateParameters;
    %% checking constraints
    cum_number_of_generation = 0;
    
     %% size2
    cum_generation_previous = cum_number_of_generation;
    hh_size = 2;
    [size2data_to_check_all,size2outcome_to_save, size2number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, cum_number_of_generation,L,ACS_count, 150000);
    
    extras_size2 = GetExtraIndividual(size2data_to_check_all,size2outcome_to_save,hh_size,ACS_count,cum_generation_previous);
    cum_number_of_generation = cum_number_of_generation + size2number_of_generation;
    
    
    %% size3
    cum_generation_previous = cum_number_of_generation;
    hh_size = 3;
    [size3data_to_check_all,size3outcome_to_save, size3number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, cum_number_of_generation,L,ACS_count,500000);
    
    extras_size3 = GetExtraIndividual(size3data_to_check_all,size3outcome_to_save,hh_size,ACS_count,cum_generation_previous);
    cum_number_of_generation = cum_number_of_generation + size3number_of_generation;
    
    %% size4
    cum_generation_previous = cum_number_of_generation;
    hh_size = 4;
    
    [size4data_to_check_all,size4outcome_to_save, size4number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, cum_number_of_generation,L,ACS_count,5000000);
    
    extras_size4 = GetExtraIndividual(size4data_to_check_all,size4outcome_to_save,hh_size,ACS_count,cum_generation_previous);
    cum_number_of_generation = cum_number_of_generation + size3number_of_generation;
    %%   
    postsave;
    
end

toc    
