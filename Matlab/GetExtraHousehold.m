function extra = GetExtraHousehold(data_to_check_all,outcome_to_save,hh_size,ACS_count)
    %good household
    extras_old2 = data_to_check_all(outcome_to_save==0,1);
    
    %bad household
    syndatao_old = data_to_check_all(size2outcome_to_save==1,1);
    exceeding_index = syndatao_old(ACS_count(hh_size-1),1);
    while (isempty(find(extras_old2(:,1)==exceeding_index))==1)
    	exceeding_index = exceeding_index - 1;
    end
    exceeding_index = exceeding_index-10000;
    extras_old1 = data_to_check_all(1:exceeding_index,:);
    
    extras_old = extras_old1(find(outcome_to_save(1:exceeding_index)==0),:);
    size_extras_old = size(extras_old),1;

    extra = zeros(hh_size*size_extras_old,8+hh_size);
    
    hh_index = hh_size * 8 + 1;
    for s = 1:size_extras_old
        ibase = 2*(s-1);
        for hm = 1: hh_size
            extra(ibase+hm,1:8) = extras_old(s,(hm-1)*8 +(1:8));
            extra(ibase+hm,9) = extras_old(s,hh_index);
            extra(ibase+hm,10)= extras_old(s,hh_index + hm);
        end
    end
end

