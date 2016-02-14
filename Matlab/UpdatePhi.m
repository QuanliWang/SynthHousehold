function [phi,phicountcluster,kcount] = UpdatePhi(data_full_all,K,L,p,d,maxd,...
    z_HH_Individuals_all,z_Individual_all,z_HH_all)

    phi=zeros(maxd,p,K*L);      % cell probabilities
    data = data_full_all(:,3:7);
       
    groupIndex = L*(z_HH_Individuals_all-1)+z_Individual_all;
    for j = 1:p
        phicount = groupcount(groupIndex,data(:,j),K*L, d(j));
        for k = 1:K
            for l = 1:L
                group = L*(k-1)+l;
                phi1 = gamrnd(phicount(group,:)+1,1); % variable j, z_HH k, z_member L
                phi(1:d(j),j,group) = phi1 / sum(phi1);
            end
        end        
    end
    
    phicountcluster = zeros(K,L);
    for k = 1:K
        zh1 = (z_HH_Individuals_all==k);
        for l = 1:L
            zh2 = (z_Individual_all==l);
            phicountcluster(k,l) = sum(zh1&zh2);
        end
    end

    levelk = 1:K;
    kcount = sum(hist(z_HH_all,levelk),1);

    disp('phi updated');
end

