function data_to_check = GenerateData2Check(hh_size,lambda, w, ...
    phi,pi, d, p, current_number_of_generation, howmany)
    
    L = size(w,2);
    data_to_check = zeros(howmany,8*hh_size+1+hh_size);
    
    %no need to normalize it as it is cared in randomsample_new function
    hhindexh = randomsample_new(pi.*lambda{2}(:,hh_size-1),rand(howmany,1));
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
    
    phi = reshape(phi,[],p,size(phi,2)); %DIM(phi) = (maxd,p,K*L) or (maxd*p, K*L)
    for hh = 1:hh_size
        group_index = (hhindexh-1)*L+data_to_check(:,hh_size * 8 +1 + hh);
        syn = sampleindividuals(phi,d,group_index,p,rand(length(group_index),p)); 
        
        data_to_check(:,(hh-1) * 8 + 1) = (1:howmany) + current_number_of_generation*howmany;
        data_to_check(:,(hh-1) * 8 + 2) = hh;  
        data_to_check(:, (hh-1) * 8 + (3:7)) = syn;
        
    end
end

