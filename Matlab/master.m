function dummy = master()
addpath dist

n = 10000;
[SS,origdata] = PrepareData(n);

%%
% variables: HHindex 1, pernum 2, sex 3, race 4, ethn 5, age 6, relate 7(ind level) ownership 8 (binary, HH level)
[n_s,p] = size(origdata);
p = p-4;

HHrowIndex = [1 cumsum(SS)+1]; 
HHdataorig = origdata(HHrowIndex(1: (end -1)),8:9);
HHdata1 = HHdataorig(:,1);
HHdata2 = HHdataorig(:,2)-1;
HHserial = origdata(:,1);
clear HHrowIndex

% variables: HHindex 1, pernum 2, sex 3, race 4, ethn 5, age 6, relate
% 7(ind level) ownership 8 (binary, HH level) household size 

d = [2,9,5,94,12]; %number of levels for each variable (5 variables all together)
data = origdata(:,3:7);
dataT = origdata(:,3:7)';
maxd = max(d);

% what are those?
dHH1 = 2;
dHH2 = 3;

%%
%%%%---global parameters---%%

%nrun=10000; burn=8000; thin=10; %50;
nrun=25; burn=20; thin=1;
eff.sam=(nrun-burn)/thin;

K=40;                     
L=15;
aa=0.25; ab=0.25;                  % gamma hyperparameters for alpha
ba=0.25; bb=0.25;

%%%%---define output files---%%
alphaout = zeros(nrun,1); 
betaout = zeros(nrun,1);
piout = zeros(eff.sam,K); 
wout = zeros(eff.sam,K,L);
nout = zeros(nrun,1);
n_sout = zeros(nrun,1);

size2extrasize = zeros(nrun,1);
size3extrasize = zeros(nrun,1);
size4extrasize = zeros(nrun,1);

z_HH_save = zeros(nrun,n_s);
z_member_save = zeros(nrun,n_s);
elapsed_time = zeros(1,nrun);

newphiout = zeros(eff.sam,maxd*p,K*L);
lambda1out = zeros(eff.sam,K,dHH1);
lambda2out = zeros(eff.sam,K,dHH2);

%%%%---initial values---%%
alpha=1;        % hyperparameters for stick-breaking weights
beta=1;

%initialize phi, rewrite as a 3D matrix to save space
% from 94*5 * K * L to sum(d) * (K*L)
% this will save 74% of space for this particular one
phi=zeros(K*L,p,maxd);      % cell probabilities
for j=1:p
    for c=1:d(j)
        phi(1,j,c)=sum(data(:,j)==c)/n_s;
    end
end
for count = 2:K*L
    phi(count,:,:) = phi(1,:,:);
end

%initialize lambda
lambda1 = zeros(K,dHH1);
for c=1:dHH1
    lambda1(:,c) = sum(HHdata1==c)/n;
end

lambda2 = zeros(K,dHH2);
for c=1:dHH2
    lambda2(:,c) = sum(HHdata2==c)/n;
end

u=[betarnd(1,alpha,[K-1,1]);1];
pi=zeros(K,1);
pi(1:K)=u.*cumprod([1;1-u(1:K-1)]);

v=[betarnd(1,beta,[K,L-1]),ones(K,1)];
w=zeros(K,L);
v1 = (v(1,:))';
w(1,1:L)=v1.*cumprod([1;1-v1(1:L-1)]);
for k=1:K
    v1 = (v(k,:))';
    w(k,1:L)=v1.*cumprod([1;1-v1(1:L-1)]);
end

ACS_count = zeros(3,1);
ACS_count(1) = sum(SS==2);
ACS_count(2) = sum(SS==3);
ACS_count(3) = sum(SS==4);

%profile on
%%%%---Blocked sampler---%%
for i = 1:nrun    
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
            
       if (i==1)
           data_full = origdata;
           HHdata1 = HHdataorig(:,1);
           HHdata2 = HHdataorig(:,2)-1;
       else
           trueextras_size2 = extras_size2;
           trueextras_size3 = extras_size3;
           trueextras_size4 = extras_size4;
           trueextras = [trueextras_size2;trueextras_size3;trueextras_size4];
           
           trueextras_HHindex = unique(trueextras(:,1));
           trueextras_SS = hist(trueextras(:,1),trueextras_HHindex); 
           trueextras_ownership = zeros(size_extras_size2_old+size_extras_size3_old+size_extras_size4_old,1);
           z_HH_new = zeros(size_extras_size2_old+size_extras_size3_old+size_extras_size4_old,1);
           trueextras_HHdatafull = trueextras(:,8);
           trueextras_HH = trueextras(:,9);
           
           cumSS = cumsum(trueextras_SS);
           trueextras(1:cumSS(1),1) = 10001;
           for s=2:size_extras_size2_old+size_extras_size3_old+size_extras_size4_old
                trueextras(cumSS(s-1)+1:cumSS(s),1) = s+10000;
                trueextras_ownership(s) = trueextras_HHdatafull(cumSS(s));
                z_HH_new(s) = trueextras_HH(cumSS(s));
           end
           data_full = [origdata(:,1:8);trueextras(:,1:8)];
           HHdata1 = [HHdataorig(:,1);trueextras_ownership];
           HHdata2 = [HHdataorig(:,2)-1;ones(size_extras_size2_old,1);2*ones(size_extras_size3_old,1);3*ones(size_extras_size4_old,1)];
       end
       
 
       data = data_full(:,3:7);
       
       HHindex_new = unique(data_full(:,1)); 
       SS_new = hist(data_full(:,1),HHindex_new); 
       
       SSsize_new = size(SS_new);
       n_new = SSsize_new(2);
       n_s_new0 = size(data);
       n_s_new = n_s_new0(1);

       z_HH= zeros(n,1);                   %cluster indicators for household; K
       z_member = zeros(n_s,1);             %cluster indicators for household member; L
       
       z_HH_prob = samplezHHwithHHnewv1_2HHvar(newphi,dataT,w,pi,SS',K,L,p,maxd,n,HHdata1(1:n),lambda1,HHdata2(1:n),lambda2);
       for h = 1:n
           zupdateprob_h = z_HH_prob(K*(h-1)+1:K*h);
           z_HH(h) = randomsample(zupdateprob_h,rand);
       end
       
       z_HH_all = zeros(n_s,1);
       cumsumSS = [0 cumsum(SS)];
       for h = 1:n
           z_HH_all(cumsumSS(h)+1:cumsumSS(h+1)) = z_HH(h);
       end
       
       if (i==1)
           z_HH_all_every = z_HH_all;
           z_HH_every = z_HH;
       else
           z_HH_all_every = [z_HH_all;trueextras(:,9)];
           z_HH_every = [z_HH;z_HH_new];
       end
       
       disp('zHH updated');
       %% update zmember
       %%% to work in C++
       
       z_member_prob = samplezmemberv1(newphi,dataT,w,K,L,p,maxd,n_s,z_HH,HHserial);

       for m = 1:n_s
           zupdateprob_m = z_member_prob(L*(m-1)+1:L*m);
           z_member(m) = randomsample(zupdateprob_m,rand);
       end

       if (i==1)
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
                value = (HHdata1==c);
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
                value = (HHdata2==c);
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
   
    %% checking constraints
    cum_number_of_generation = 0;
     %% size2
    hh_size = 2;
    size2data_to_check_all = zeros(150000,8*hh_size+1+hh_size);
    size2outcome_to_save = zeros(150000,1);
    
    [size2data_to_check_all,size2outcome_to_save, size2number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, cum_number_of_generation,L,size2outcome_to_save,size2data_to_check_all,ACS_count);


    syndatao_size2_old = size2data_to_check_all(find(size2outcome_to_save==1),:);
    extras_size2_old2 = size2data_to_check_all(find(size2outcome_to_save==0),:);
   
    syndatao_size2 = syndatao_size2_old(1:ACS_count(hh_size-1),:);
    size2exceeding_index = syndatao_size2(ACS_count(hh_size-1),1);
    
    while (isempty(find(extras_size2_old2(:,1)==size2exceeding_index))==1)
    	size2exceeding_index = size2exceeding_index - 1;
    end
    size2exceeding_index = size2exceeding_index-10000;
    extras_size2_old1 = size2data_to_check_all(1:size2exceeding_index,:);
    extras_size2_old = extras_size2_old1(find(size2outcome_to_save(1:size2exceeding_index)==0),:);
    
    size_extras_size2_old_1 = size(extras_size2_old);
    size_extras_size2_old = size_extras_size2_old_1(1);

    extras_size2 = zeros(2*size_extras_size2_old,8+2);
    for s = 1:size_extras_size2_old
        extras_size2(2*(s-1)+1,1:8) = extras_size2_old(s,1:8);
        extras_size2(2*(s-1)+2,1:8) = extras_size2_old(s,9:16);
        extras_size2((2*(s-1)+1):2*s,9) = extras_size2_old(s,17);
        extras_size2(2*(s-1)+1,10) = extras_size2_old(s,18);
        extras_size2(2*(s-1)+2,10) = extras_size2_old(s,19);
    end
    
    cum_number_of_generation = cum_number_of_generation + size2number_of_generation;
    
    %% size3
    hh_size = 3;
    size3data_to_check_all = zeros(500000,8*hh_size+1+hh_size);
    size3outcome_to_save = zeros(500000,1);
    [size3data_to_check_all,size3outcome_to_save, size3number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, cum_number_of_generation,L,size3outcome_to_save,size3data_to_check_all,ACS_count);

    syndatao_size3_old = size3data_to_check_all(find(size3outcome_to_save==1),:);
    extras_size3_old2 = size3data_to_check_all(find(size3outcome_to_save==0),:);
   
    syndatao_size3 = syndatao_size3_old(1:ACS_count(hh_size-1),:);
    size3exceeding_index = syndatao_size3(ACS_count(hh_size-1),1);

    while (isempty(find(extras_size3_old2(:,1)==size3exceeding_index))==1)
    	size3exceeding_index = size3exceeding_index - 1;
    end
    size3exceeding_index = size3exceeding_index-(1+size2number_of_generation)*10000;
    extras_size3_old1 = size3data_to_check_all(1:size3exceeding_index,:);
    extras_size3_old = extras_size3_old1(find(size3outcome_to_save(1:size3exceeding_index)==0),:);
    
    size_extras_size3_old_1 = size(extras_size3_old);
    size_extras_size3_old = size_extras_size3_old_1(1);

    extras_size3 = zeros(3*size_extras_size3_old,8+2);
    for s = 1:size_extras_size3_old
        extras_size3(3*(s-1)+1,1:8) = extras_size3_old(s,1:8);
        extras_size3(3*(s-1)+2,1:8) = extras_size3_old(s,9:16);
        extras_size3(3*(s-1)+3,1:8) = extras_size3_old(s,17:24);
        extras_size3((3*(s-1)+1):3*s,9) = extras_size3_old(s,25);
        extras_size3(3*(s-1)+1,10) = extras_size3_old(s,26);
        extras_size3(3*(s-1)+2,10) = extras_size3_old(s,27);
        extras_size3(3*(s-1)+3,10) = extras_size3_old(s,28);
    end
    cum_number_of_generation = cum_number_of_generation + size3number_of_generation;
    
    %% size4
    hh_size = 4;
    size4data_to_check_all = zeros(5000000,8*hh_size+1+hh_size);
    size4outcome_to_save = zeros(5000000,1);
    
    [size4data_to_check_all,size4outcome_to_save, size4number_of_generation] = ...
    GenerateData2Check_all(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, cum_number_of_generation,L,size4outcome_to_save,size4data_to_check_all,ACS_count);
    
    syndatao_size4_old = size4data_to_check_all(find(size4outcome_to_save==1),:);
    extras_size4_old2 = size4data_to_check_all(find(size4outcome_to_save==0),:);
   
    syndatao_size4 = syndatao_size4_old(1:ACS_count(hh_size-1),:);
    size4exceeding_index = syndatao_size4(ACS_count(hh_size-1),1);
    
    while (isempty(find(extras_size4_old2(:,1)==size4exceeding_index))==1)
    	size4exceeding_index = size4exceeding_index - 1;
    end
    size4exceeding_index = size4exceeding_index-(1+size2number_of_generation+size3number_of_generation)*10000;
    extras_size4_old1 = size4data_to_check_all(1:size4exceeding_index,:);
    extras_size4_old = extras_size4_old1(find(size4outcome_to_save(1:size4exceeding_index)==0),:);
    
    size_extras_size4_old_1 = size(extras_size4_old);
    size_extras_size4_old = size_extras_size4_old_1(1);

    extras_size4 = zeros(4*size_extras_size4_old,8+2);
    for s = 1:size_extras_size4_old
        extras_size4(4*(s-1)+1,1:8) = extras_size4_old(s,1:8);
        extras_size4(4*(s-1)+2,1:8) = extras_size4_old(s,9:16);
        extras_size4(4*(s-1)+3,1:8) = extras_size4_old(s,17:24);
        extras_size4(4*(s-1)+4,1:8) = extras_size4_old(s,25:32);
        extras_size4((4*(s-1)+1):4*s,9) = extras_size4_old(s,33);
        extras_size4(4*(s-1)+1,10) = extras_size4_old(s,34);
        extras_size4(4*(s-1)+2,10) = extras_size4_old(s,35);
        extras_size4(4*(s-1)+3,10) = extras_size4_old(s,36);
        extras_size4(4*(s-1)+4,10) = extras_size4_old(s,37);
    end
    
    
    %%    
    disp('synthesis updated');
    
       %syndata{i} = syndatao;
       if (mod(i,thin) == 0 && i > burn) 
           piout((i-burn)/thin,:) = pi';
           wout((i-burn)/thin,:,:) = w;
           newphiout((i-burn)/thin,:,:) = newphi;
           lambda1out((i-burn)/thin,:,:) = lambda1;
           lambda2out((i-burn)/thin,:,:) = lambda2;
       end
       i
       n_new
       toc
       elapsed_time(i) = toc;
       n_sout(i) = n_s_new;
       nout(i) = n_new;
       size2extrasize(i) = size_extras_size2_old;
       size3extrasize(i) = size_extras_size3_old;
       size4extrasize(i) = size_extras_size4_old;
       z_HH_save(i,1:n_s) = z_HH_all;
       z_member_save(i,1:n_s) = z_member;
       alphaout(i) = alpha;
       betaout(i) = beta;
       disp('saving done'); 
end

toc    
