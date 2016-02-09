%retuen extra impossible households
function [Individuals,HouseHolds] = GetExtraIndividualAndHousehold(...
    data_to_check_all,outcome_to_save,hh_size,ACS_count)
    %impossible household
    extras_old2 = data_to_check_all(outcome_to_save==0,:);

    %possible household
    possiblehh = data_to_check_all(outcome_to_save==1,1);
    
    %find the number ACS_count(hh_size-1)'th possible household 
    exceeding_index = possiblehh(ACS_count(hh_size-1),1);
    
    %Get all the impossible households we need
    HouseHolds = extras_old2(extras_old2(:,1) < exceeding_index ...
        & extras_old2(:,1) >0,:);
     
    extra_hh = size(HouseHolds,1);
    Individuals = zeros(hh_size*extra_hh,8+2);
    hh_index = hh_size * 8 + 1;
    temp = 1:8;
    for s = 1:extra_hh
        ibase = hh_size*(s-1);
        for hm = 1: hh_size
            Individuals(ibase+hm,1:8) = HouseHolds(s, (hm-1)*8 + temp);
            Individuals(ibase+hm,9) = HouseHolds(s,hh_index);
            Individuals(ibase+hm,10)= HouseHolds(s,hh_index + hm);
        end
    end 
end