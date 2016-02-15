function [lambda1,lambda2] = UpdateLambda(dHH1,dHH2,K,z_HH_all,HHdata1_all,HHdata2_all)
    lambda1 = zeros(K,dHH1);
    lambda1count = zeros(K,dHH1);
    for k = 1:K
        lam1 = (z_HH_all==k); 
        for c = 1:dHH1
            value = (HHdata1_all==c);
            lambda1count(k,c) = lambda1count(k,c)+sum(lam1&value);
        end
        lam11 = gamrnd(lambda1count(k,1:dHH1)+1,1); % variable j, z_HH k, z_member L
        lam111 = bsxfun(@times,lam11(1:dHH1),1./sum(lam11(1:dHH1)));
        lambda1(k,1:dHH1) = lam111;
    end

    lambda2 = zeros(K,dHH2);
    lambda2count = zeros(K,dHH2);
    for k = 1:K
        lam1 = (z_HH_all==k); 
        for c = 1:dHH2
            value = (HHdata2_all==c);
            lambda2count(k,c) = lambda2count(k,c)+sum(lam1&value);
        end
        lam11 = gamrnd(lambda2count(k,1:dHH2)+1,1); % variable j, z_HH k, z_member L
        lam111 = bsxfun(@times,lam11(1:dHH2),1./sum(lam11(1:dHH2)));
        lambda2(k,1:dHH2) = lam111;
    end
end

