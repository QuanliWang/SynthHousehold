%retuen extra impossible households
function [Individuals,HouseHolds] = GetExtraIndividualAndHousehold(...
    data_to_check_all,outcome_to_save,hh_size,ACS_count)
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