function syn = GenerateIndividualVariables(phi,p,d,group_index)
    syn = zeros(length(group_index),p);
    levelGroup = unique(group_index);
    groupCount = hist(group_index,levelGroup);
    for j = 1:p    
        phi_j = squeeze(phi(1:d(j),j,:));

        for mg = 1: length(levelGroup)
            syn(group_index == levelGroup(mg),j) = ...
                randomsample_new(phi_j(:,levelGroup(mg)),rand(groupCount(mg),1));
        end
    end
end

