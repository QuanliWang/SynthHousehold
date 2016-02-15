function [lambda1,lambda2] = UpdateLambda(dHH1,dHH2,K,z_HH_all,HHdata1_all,HHdata2_all)

    lambda1count = groupcount(z_HH_all,HHdata1_all,K, dHH1);    
    lambda1 = zeros(K,dHH1);
    for k = 1:K
        lam11 = gamrnd(lambda1count(k,1:dHH1)+1,1); % variable j, z_HH k, z_member L
        lambda1(k,:) = lam11/sum(lam11);
    end
    
    lambda2 = zeros(K,dHH2);
    lambda2count = groupcount(z_HH_all,HHdata2_all,K, dHH2); 
    for k = 1:K
        lam11 = gamrnd(lambda2count(k,1:dHH2)+1,1); % variable j, z_HH k, z_member L
        lambda2(k,:) = lam11/sum(lam11);
    end

end

