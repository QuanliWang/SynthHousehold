function [phi,phicountcluster,kcount,phicount] = UpdatePhi(data_full_all,K,L,p,d,maxd,...
    z_HH_Individuals_all,z_Individual_all,z_HH_all)
    %%      
    % -- update phi-- %
    %tic;
    phi=zeros(K*L,p,maxd);      % cell probabilities
    data = data_full_all(:,3:7);
    phicount = zeros(K*L,p,maxd);
    for m = 1:size(data,1)
        mdata = data(m,:);
        for j = 1:p
            mk = z_HH_Individuals_all(m); 
            lk = z_Individual_all(m);
            phicount(L*(mk-1)+lk,j,mdata(j)) = phicount(L*(mk-1)+lk,j,mdata(j))+1;
        end
    end

    for j = 1:p
        phicountp = reshape(phicount(:,j,:),K*L,maxd);
        for k = 1:K
            for l = 1:L
                phi1 = gamrnd(phicountp(L*(k-1)+l,1:d(j))+1,1); % variable j, z_HH k, z_member L
                phi11 = bsxfun(@times,phi1(1:d(j)),1./sum(phi1(1:d(j))));
                phi(L*(k-1)+l,j,1:d(j)) = phi11;
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

