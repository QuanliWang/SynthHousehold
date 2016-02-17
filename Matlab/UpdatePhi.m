function [phi,phicountcluster,kcount] = UpdatePhi(K,L,p,d,maxd,...
    IndividualData_all, z_Individual_all,z_HH_all,z_HH_Individuals_all)

    phi=zeros(maxd,p,K*L);      % cell probabilities
    data = IndividualData_all(:,3:7);
       
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
    phi = reshape(phi,maxd*p, K * L); %reshape to a 2D matrix 
    phicountcluster = groupcount(z_HH_Individuals_all,z_Individual_all,K,L);
    
    levelk = 1:K;
    kcount = sum(hist(z_HH_all,levelk),1);
end

