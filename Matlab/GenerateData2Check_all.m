function [data_to_check_all,outcome_to_save, number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, cum_number_of_generation,L,outcome_to_save,data_to_check_all,ACS_count)
    number_of_generation = 0;
    while (sum(outcome_to_save)<=ACS_count(hh_size-1))
        number_of_generation = number_of_generation + 1;
        data_to_check = GenerateData2Check(hh_size,lambda1, lambda2, w, ...
            phi,pi, d, p, number_of_generation,cum_number_of_generation, L);
        
        % transpose the data to check
        data_to_checkT = data_to_check(:,1:(hh_size*8))';   
        if hh_size == 2 
            outcome = checkingconstraints_size2_sorted(data_to_checkT,10000);
        elseif hh_size == 3
            outcome = checkingconstraints_size3_sorted(data_to_checkT,10000);
        elseif hh_size == 4
            outcome1 = checkingconstraints_size4_sorted_part1(data_to_checkT,10000);
            outcome2 = checkingconstraints_size4_sorted_part2(data_to_checkT,10000);
            outcome3 = checkingconstraints_size4_sorted_part3(data_to_checkT,10000);
            outcome = outcome1+outcome2+outcome3;
        else 
            error( ['Household size ' num2str(hh_size) ' not considered']);
        end
            
        startindex1 = min(find(data_to_check_all(:,1)==0));
        startindex2 = (number_of_generation-1)*10000 + 1;
        data_to_check_all(startindex1:startindex1+10000-1,:) = data_to_check;
        outcome_to_save(startindex2:(startindex2+10000-1)) = outcome;
    end
end

