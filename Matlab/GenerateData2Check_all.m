function [Individuals,HouseHolds, number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda, w, ...
    phi,pi, d, p, cum_number_of_generation,L,ACS_count,newxtras)
    
    %Individuals_all = [];
    %HouseHolds_all = [];
    data_to_check_all = zeros(newxtras,8*hh_size+1+hh_size);
    outcome_to_save = zeros(newxtras,1);
    number_of_generation = 0;
    while (sum(outcome_to_save)<=ACS_count(hh_size-1))
        number_of_generation = number_of_generation + 1;
        
        %generate a batch of 10K household
        data_to_check = GenerateData2Check(hh_size,lambda, w, ...
            phi,pi, d, p, number_of_generation+cum_number_of_generation, L);
        
        % transpose the data to check
        data_to_checkT = data_to_check(:,1:(hh_size*8))';   
        
        outcome = checkingconstraints(data_to_checkT,10000,hh_size);
        
        startindex = (number_of_generation-1)*10000 + 1;
        endindex = (startindex+10000-1);
        data_to_check_all(startindex:endindex,:) = data_to_check;
        outcome_to_save(startindex:endindex) = outcome;
    end
    
    number_of_generation = number_of_generation + cum_number_of_generation;
    
    %impossible household
    extras_old2 = data_to_check_all(outcome_to_save==0 &...
        data_to_check_all(:,1) > 0,:);

    %possible household
    possiblehh = data_to_check_all(outcome_to_save==1,1);
    
    %find the number ACS_count(hh_size-1)'th possible household 
    exceeding_index = possiblehh(ACS_count(hh_size-1),1);
    
    %Get all the impossible households we need
    HouseHolds = extras_old2(extras_old2(:,1) < exceeding_index ...
        & extras_old2(:,1) >0,:);
     
    t1 = HouseHolds(:,1:hh_size * 8)';
    t2 = HouseHolds(:,hh_size * 8 +1)';
    t3 = HouseHolds(:,hh_size * 8 +1 + (1:hh_size))';
    
    Individuals = [reshape(t1,8,[])' reshape(repmat(t2,hh_size,1),[],1)...
        reshape(t3,[],1)];
    
end
