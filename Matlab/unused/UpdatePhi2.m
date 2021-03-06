function phi = UpdatePhi2(K,L,p,d,maxd, IndividualData, z_Individual,picount)
    phi=zeros(maxd,p,K*L);      % cell probabilities
    data = IndividualData(:,3:7);
    groupIndex = L*(z_Individual(:,1)-1)+z_Individual(:,2);
    for j = 1:p
        phicount_all = picount(:,1:d(j),j) + groupcount(groupIndex,data(:,j),K*L, d(j));
        for k = 1:K
            for l = 1:L
                group = L*(k-1)+l;
                phi1 = gamrnd(phicount_all(group,:)+1,1); % variable j, z_HH k, z_member L
                phi(1:d(j),j,group) = phi1 / sum(phi1);
            end
        end        
    end
    phi = reshape(phi,maxd*p, K * L); %reshape to a 2D matrix 
end

