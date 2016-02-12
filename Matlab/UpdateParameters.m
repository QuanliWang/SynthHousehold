tic;   
        %% update zHH
       %%% to work in C++, consider make phi into a 2-d matrix maxd*p-by-K*L
       newphi = zeros(maxd*p,K*L);
       for e = 1:K*L
           newphimatrix = reshape(phi(e,:,:),p,maxd);
           for j = 1:p
               newphi(maxd*(j-1)+1:maxd*j,e) = newphimatrix(j,:);
           end
       end
            
       data = data_full_all(:,3:7);
       
       HHindex_new = unique(data_full_all(:,1)); 
       SS_new = hist(data_full_all(:,1),HHindex_new); 
       
       SSsize_new = size(SS_new);
       n_new = size(SS_new,2);
       n_s_new0 = size(data);
       n_s_new = n_s_new0(1);

       z_HH= zeros(n,1);                   %cluster indicators for household; K
       z_member = zeros(n_s,1);             %cluster indicators for household member; L
       
       z_HH_prob = samplezHHwithHHnewv1_2HHvar(newphi,dataT,w,pi,SS',K,L,p,maxd,n,HHdata1_all(1:n),lambda1,HHdata2_all(1:n),lambda2);
       for h = 1:n
           zupdateprob_h = z_HH_prob(K*(h-1)+1:K*h);
           z_HH(h) = randomsample(zupdateprob_h,rand);
       end
       
       z_HH_all = zeros(n_s,1);
       cumsumSS = [0 cumsum(SS)];
       for h = 1:n
           z_HH_all(cumsumSS(h)+1:cumsumSS(h+1)) = z_HH(h);
       end
       
       if isempty(trueextras)
            z_HH_all_every = z_HH_all;
       else
            z_HH_all_every = [z_HH_all;trueextras(:,9)];
       end
       
       z_HH_every = [z_HH;z_HH_extra];
       disp('zHH updated');

       %% update zmember
       %%% to work in C++
       z_member_prob = samplezmemberv1(newphi,dataT,w,K,L,p,maxd,n_s,z_HH,HHserial);

       for m = 1:n_s
           zupdateprob_m = z_member_prob(L*(m-1)+1:L*m);
           z_member(m) = randomsample(zupdateprob_m,rand);
       end
       
       if isempty(trueextras)
           z_member_every = z_member;
       else
           z_member_every = [z_member;trueextras(:,10)];
       end
       
       %toc
       disp('zmember updated');
       %%      
        % -- update phi-- %
        %tic;
        phicount = zeros(K*L,p,maxd);
        for m = 1:n_s_new
            mdata = data(m,:);
            for j = 1:p
                mk = z_HH_all_every(m); 
                lk = z_member_every(m);
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
            zh1 = (z_HH_all_every==k);
            for l = 1:L
                zh2 = (z_member_every==l);
                phicountcluster(k,l) = sum(zh1&zh2);
            end
        end
     
        levelk = 1:K;
        kcount = sum(hist(z_HH_every,levelk),1);
        
        disp('phi updated');

        %%
        % -- update lambda -- %
        %tic;
        lambda1count = zeros(K,dHH1);
        for k = 1:K
            lam1 = (z_HH_every==k); 
            for c = 1:dHH1
                value = (HHdata1_all==c);
                lambda1count(k,c) = lambda1count(k,c)+sum(lam1&value);
            end
            lam11 = gamrnd(lambda1count(k,1:dHH1)+1,1); % variable j, z_HH k, z_member L
            lam111 = bsxfun(@times,lam11(1:dHH1),1./sum(lam11(1:dHH1)));
            lambda1(k,1:dHH1) = lam111;
        end
        
        lambda2count = zeros(K,dHH2);
        for k = 1:K
            lam1 = (z_HH_every==k); 
            for c = 1:dHH2
                value = (HHdata2_all==c);
                lambda2count(k,c) = lambda2count(k,c)+sum(lam1&value);
            end
            lam11 = gamrnd(lambda2count(k,1:dHH2)+1,1); % variable j, z_HH k, z_member L
            lam111 = bsxfun(@times,lam11(1:dHH2),1./sum(lam11(1:dHH2)));
            lambda2(k,1:dHH2) = lam111;
        end
        
        %toc
        disp('lambda updated');
       %%         
        % -- update pi -- % 
        %tic;
        for k = 1:K-1
            u(k) = betarnd(1 + kcount(k),alpha + sum(kcount(k+1:K)));
            if u(k)>1-1e-5
                u(k)=1-1e-5;
            end  
        end
        u(K) = 1;
        
        pi(1:K)=u.*cumprod([1;1-u(1:K-1)]);
        
        %toc
        disp('pi updated');
        % -- update w -- %
        %tic;       
        for k = 1:K
            vk = v(k,:);
            for l = 1:L-1
                vk(l) = betarnd(1 + phicountcluster(k,l), beta + sum(phicountcluster(k,l+1:L)));
                if vk(l)>1-1e-5
                    vk(l)=1-1e-5;
                end
                vk(L) = 1;
            end
            v(k,:) = vk;
            w(k,1:L) = vk.*cumprod([1,1-vk(1:L-1)]);
        end
        %toc
        disp('w updated');
        %tic;
        % -- update alpha-- %
        alpha = gamrnd(aa + K - 1, 1/(ab - sum(log(1-u(1:K-1)))));
        
        % -- update beta -- %
        vnew = v(:,1:L-1);
        beta = gamrnd(ba + K*(L-1), 1/(bb - sum(sum(log(1-vnew)))));
        %toc
        
        disp('alpha beta updated');
    %tic;
   