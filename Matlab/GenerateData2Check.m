function data_to_check = GenerateData2Check(hh_size,lambda, w, ...
    phi,pi, d, p, current_number_of_generation)
    L = size(w,2);
    data_to_check = zeros(10000,8*hh_size+1+hh_size);
    
    %no need to normalize it as it is cared in randomsample_new function
    hhindexh = randomsample_new(pi.*lambda{2}(:,hh_size-1),rand(10000,1));
    data_to_check(:,hh_size * 8 + 1) = hhindexh;
    
    %memberindexhh
    %do random samples for the same probs at the same time
    levelhh = unique(hhindexh); 
    cluster_count = hist(hhindexh,levelhh);
    for i = 1:length(levelhh)
        memberindexhh = randomsample_new(w(levelhh(i),:),rand(hh_size*cluster_count(i),1));
        data_to_check(hhindexh == levelhh(i),hh_size * 8 +1 + (1:hh_size))=...
            reshape(memberindexhh,[],hh_size);
        
        % generating household level data
        data_to_check(hhindexh == levelhh(i),(1:hh_size)*8) = ...
            repmat(randomsample_new(lambda{1}(levelhh(i),:),rand(cluster_count(i),1)), 1,hh_size);
        
    end
    
    phi = reshape(phi,[],p,size(phi,2));
    for hh = 1:hh_size
        syn = zeros(10000,p);
        group_index = (hhindexh-1)*L+data_to_check(:,hh_size * 8 +1 + hh);
        levelGroup = unique(group_index);
        groupCount = hist(group_index,levelGroup);
        for j = 1:p    
            for mg = 1: length(levelGroup)
                syn(group_index == levelGroup(mg),j) = ...
                    randomsample_new(phi(1:d(j),j,levelGroup(mg)),rand(groupCount(mg),1));
            end
        end
        
        data_to_check(:, (hh-1) * 8 + (3:7)) = syn;
        data_to_check(:,(hh-1) * 8 + 1) = (1:10000) + current_number_of_generation*10000;
        data_to_check(:,(hh-1) * 8 + 2) = hh;  
    end
end

