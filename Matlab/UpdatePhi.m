function phi = UpdatePhi(K,L,p,d,maxd, IndividualData_all, z_Individual_all)

    phi=zeros(maxd,p,K*L);      % cell probabilities
    data = IndividualData_all(:,3:7);
    
    groupIndex = L*(z_Individual_all(:,1)-1)+z_Individual_all(:,2);
    for j = 1:p
        phicount = groupcount(groupIndex,data(:,j),K*L, d(j));
        
        Phi_j = gamrnd(phicount + 1, 1);
        sumPhi_j = sum(Phi_j,2);
        for group = 1: K * L
            Phi_j(group,:) = Phi_j(group,:) / sumPhi_j(group);
        end
        phi(1:d(j),j,:) = Phi_j';
    end
    phi = reshape(phi,maxd*p, K * L); %reshape to a 2D matrix 
end

