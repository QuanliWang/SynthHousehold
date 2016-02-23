function [lambda] = UpdateLambda(dHH,K,z_HH_all,HHdata_all)
    lambda = cell(2,1);
    lambda1count = groupcount(z_HH_all,HHdata_all(1,:),K, dHH(1));    
    lambda1 = zeros(K,dHH(1));
    for k = 1:K
        lam11 = gamrnd(lambda1count(k,1:dHH(1))+1,1); % variable j, z_HH k, z_member L
        lambda1(k,:) = lam11/sum(lam11);
    end
    
    lambda2 = zeros(K,dHH(2));
    lambda2count = groupcount(z_HH_all,HHdata_all(2,:),K, dHH(2)); 
    for k = 1:K
        lam11 = gamrnd(lambda2count(k,1:dHH(2))+1,1); % variable j, z_HH k, z_member L
        lambda2(k,:) = lam11/sum(lam11);
    end
    lambda{1} = lambda1;
    lambda{2} = lambda2;
end

