function data_to_check = GenerateData2Check(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, current_number_of_generation,L)

    data_to_check = zeros(10000,8*hh_size+1+hh_size);
    pi_size = pi.*lambda2(:,hh_size-1);
    pi_size_renorm = pi_size./sum(pi_size);
    hhindexh = randomsample_new(pi_size_renorm,rand(10000,1));
    data_to_check(:,hh_size * 8 + 1) = hhindexh;
    
    %memberindexhh
    %do random samples for the same probs at the same time
    levelhh = unique(hhindexh);
    synPlusOne = zeros(10000,1);
    
    cluster_count = hist(hhindexh,levelhh);
    for i = 1:length(levelhh)
        memberindexhh = randomsample_new(w(levelhh(i),:),rand(hh_size*cluster_count(i),1));
        data_to_check(hhindexh == levelhh(i),hh_size * 8 +1 + (1:hh_size))=...
            reshape(memberindexhh,[],hh_size);
        
        % generating household level data
        lambda1mj = lambda1(levelhh(i),:);
        synPlusOne(hhindexh == levelhh(i)) = ...
            randomsample_new(lambda1mj,rand(cluster_count(i),1));
    end
    
    for hh = 1:hh_size
        syn = zeros(10000,p);
        for j = 1:p
            group_index = (hhindexh-1)*L+data_to_check(:,hh_size * 8 +1 + hh);
            levelGroup = unique(group_index);
            groupCount = hist(group_index,levelGroup);
            for mg = 1: length(levelGroup)
                phimj = phi(1:d(j),j,levelGroup(mg));
                syn(group_index == levelGroup(mg),j) = ...
                    randomsample_new(phimj,rand(groupCount(mg),1));
            end
        end
        
        data_to_check(:, (hh-1) * 8 + (3:7)) = syn; 
        data_to_check(:,hh * 8) = synPlusOne;
        data_to_check(:,(hh-1) * 8 + 1) = (1:10000) + current_number_of_generation*10000;
        data_to_check(:,(hh-1) * 8 + 2) = hh;  
    end
end

