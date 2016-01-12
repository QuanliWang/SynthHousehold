function dummy = master()
% 
% mex samplezHHwithHHnewv1_2HHvar.cpp
% mex samplezmemberv1.cpp
% mex randomsample.cpp
% mex checkingconstraints_size2_sorted.cpp
% mex checkingconstraints_size3_sorted.cpp
% mex checkingconstraints_size4_sorted_part1.cpp
% mex checkingconstraints_size4_sorted_part2.cpp
% mex checkingconstraints_size4_sorted_part3.cpp

size234alldata = load('ACShouseholddata_size234all.mat');
newdata = size234alldata.newdata;

HHindex = unique(newdata(:,1)); 
SSfull = hist(newdata(:,1),HHindex); 

% for i=1:length(SSfull)
%     newdata((sum(SSfull(1:i))-SSfull(i)+1):(sum(SSfull(1:i))),1) = i;
%     i;
% end

randindex = randsample(length(SSfull),10000);
SS = SSfull(randindex);

totalcount = 10000;
origdata00 = zeros(1,8);
for i=1:totalcount
    index_i = find(newdata(:,1)==randindex(i));
    origdata00 = [origdata00;newdata(index_i,:)];
    i;
end

origdata0 = [origdata00(2:length(origdata00),:),zeros(length(origdata00)-1,1)];

for i=1:totalcount
    origdata0((sum(SS(1:i))-SS(i)+1):(sum(SS(1:i))),1) = i;
    origdata0((sum(SS(1:i))-SS(i)+1):(sum(SS(1:i))),9) = SS(i);
end
    
%ACSdata= load('ACShouseholddata_size234_10000HHs.mat');  

ACSdatasample = origdata0;
sizeACSdata1 = size(ACSdatasample);
sizeACSdata = sizeACSdata1(1);

ACSsize2_count = sum(SS==2);
ACSsize3_count = sum(SS==3);
ACSsize4_count = sum(SS==4);
%{
ACSdatasample_try = ACSdatasample;
[Y,I]=sort(ACSdatasample_try(:,1));
ACSdatasample=ACSdatasample_try(I,:); %use the column indices from sort() to sort all columns of A.
%}
%%
HHindex = unique(ACSdatasample(:,1)); % serial from 1 to n=23331, total 10000 households
SS = hist(ACSdatasample(:,1),HHindex); % sizes of each HH

% SizeSerial = [HHindex,SS];

for i=1:10000
    ACSdatasample(sum(SS(1:(i-1)))+1:sum(SS(1:i)),1) = i;
end

% variables: HHindex 1, pernum 2, sex 3, race 4, ethn 5, age 6, relate 7(ind level) ownership 8 (binary, HH level)
[n_s,p] = size(ACSdatasample);
nn = size(SS);
n = nn(2);
p = p-4;
ACS = ACSdatasample;

data = ACS(:,3:7);
HHdatafull = ACS(:,8:9);
HHdataorig = zeros(n,2);
for i=1:n
    HHdataorig(i,:) = HHdatafull(sum(SS(1:i-1))+1,:);
end

HHdata1 = HHdataorig(:,1);
HHdata2 = HHdataorig(:,2)-1;
origdata = ACS;
first2data = ACS(:,1:2);

HHserial = ACS(:,1);

% variables: HHindex 1, pernum 2, sex 3, race 4, ethn 5, age 6, relate
% 7(ind level) ownership 8 (binary, HH level) household size 
d = [2,9,5,94,12];
dHH1 = 2;
dHH2 = 3;
maxd = max(d);

%%
%%%%---global parameters---%%

nrun=10000; burn=8000; thin=10; %50; 
eff.sam=(nrun-burn)/thin;

K=40;%5 %50;                      
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
%kout = zeros(eff.sam,1); 
%lout = zeros(eff.sam,1);
z_HH_save = zeros(nrun,n_s);
z_member_save = zeros(nrun,n_s);
elapsed_time = zeros(1,nrun);
newphiout = zeros(eff.sam,maxd*p,K*L);
lambda1out = zeros(eff.sam,K,dHH1);
lambda2out = zeros(eff.sam,K,dHH2);
%%%%---initial values---%%

alpha=1;        % hyperparameters for stick-breaking weights
beta=1;

phi=zeros(K*L,p,maxd);      % cell probabilities

for l=1:L
    for j=1:p
        for c=1:d(j)
        phi(l,j,c)=sum(data(:,j)==c)/n_s;
        end
    end
end

for k=1:(K-1)
    phi((k*L+1):(k*L+L),:,:) = phi(1:L,:,:);
end

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
           sizeextras1 = size(trueextras);
           sizeextras = sizeextras1(1);
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
       
       HHserial_new = data_full(:,1);
       
       data = data_full(:,3:7);
       
       %size_data1 = size(data);
       %size_data = size_data1(1);
       % also need to transpose data
       dataT = origdata(:,3:7)';
       
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
     %% size2
    size2data_to_check_all = zeros(150000,8*2+1+2);
    size2outcome_to_save = zeros(150000,1);
    size2number_of_generation = 0;
    while (sum(size2outcome_to_save)<=ACSsize2_count)
    size2number_of_generation = size2number_of_generation + 1;
    %tic;
    size2data_to_check = zeros(10000,8*2+1+2);
    for m=1:10000
        pi_size2 = pi.*lambda2(:,1);
        pi_size2_renorm = pi_size2./sum(pi_size2);
        hhindexh = randomsample(pi_size2_renorm,rand);
        size2data_to_check(m,17) = hhindexh;
        syn = zeros(2,p+1);
        for hh = 1:2
            memberindexhh = randomsample(w(hhindexh,:),rand);
            size2data_to_check(m,17+hh) = memberindexhh;
            % generating individual level data
            for j = 1:p
                phimj = reshape(phi((hhindexh-1)*L+memberindexhh,j,1:d(j)),1,d(j));
                syn(hh,j) = randomsample(phimj,rand);
            end
        end
        % generating household level data
        lambda1mj = lambda1(hhindexh,:);
        syn(:,p+1) = randomsample(lambda1mj,rand);
        syn_sorted = sortrows(syn,5);
        size2data_to_check(m,3:8) = syn_sorted(1,:);
        size2data_to_check(m,11:16) = syn_sorted(2,:);
        size2data_to_check(m,1) = m + (size2number_of_generation)*10000;
        size2data_to_check(m,9) = m + (size2number_of_generation)*10000;
        size2data_to_check(m,2) = 1;
        size2data_to_check(m,10) = 2;
        size2data_to_check(m,17) = hhindexh;
    end
    
    % transpose the data to check
    size2data_to_checkT = size2data_to_check(:,1:16)';
    %toc   
    %mex checkingconstraints_new.cpp       
    size2outcome = checkingconstraints_size2_sorted(size2data_to_checkT,10000);
    size2startindex1 = min(find(size2data_to_check_all(:,1)==0));
    % startindex2 = min(find(outcome_to_save==0));
    size2startindex2 = (size2number_of_generation-1)*10000 + 1;
    size2data_to_check_all(size2startindex1:size2startindex1+10000-1,:) = size2data_to_check;
    size2outcome_to_save(size2startindex2:(size2startindex2+10000-1)) = size2outcome;
    %sum(outcome)
    end
    
    syndatao_size2_old = size2data_to_check_all(find(size2outcome_to_save==1),:);
    extras_size2_old2 = size2data_to_check_all(find(size2outcome_to_save==0),:);
   
    syndatao_size2 = syndatao_size2_old(1:ACSsize2_count,:);
    size2exceeding_index = syndatao_size2(ACSsize2_count,1);
    
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
    
    %% size3
    size3data_to_check_all = zeros(500000,8*3+1+3);
    size3outcome_to_save = zeros(500000,1);
    size3number_of_generation = 0;
    while (sum(size3outcome_to_save)<=ACSsize3_count)
    size3number_of_generation = size3number_of_generation + 1;
    %tic;
    size3data_to_check = zeros(10000,8*3+1+3);
    for m=1:10000
        pi_size3 = pi.*lambda2(:,2);
        pi_size3_renorm = pi_size3./sum(pi_size3);
        hhindexh = randomsample(pi_size3_renorm,rand);
        size3data_to_check(m,25) = hhindexh;
        syn = zeros(3,p+1);
        for hh = 1:3
            memberindexhh = randomsample(w(hhindexh,:),rand);
            size3data_to_check(m,25+hh) = memberindexhh;
            % generating individual level data
            for j = 1:p
                phimj = reshape(phi((hhindexh-1)*L+memberindexhh,j,1:d(j)),1,d(j));
                syn(hh,j) = randomsample(phimj,rand);
            end
        end
        % generating household level data
        lambda1mj = lambda1(hhindexh,:);
        syn(:,p+1) = randomsample(lambda1mj,rand);
        syn_sorted = sortrows(syn,5);
        size3data_to_check(m,3:8) = syn_sorted(1,:);
        size3data_to_check(m,11:16) = syn_sorted(2,:);
        size3data_to_check(m,19:24) = syn_sorted(3,:);
        size3data_to_check(m,1) = m + (size2number_of_generation+size3number_of_generation)*10000;
        size3data_to_check(m,9) = m + (size2number_of_generation+size3number_of_generation)*10000;
        size3data_to_check(m,17) = m + (size2number_of_generation+size3number_of_generation)*10000;
        size3data_to_check(m,2) = 1;
        size3data_to_check(m,10) = 2;
        size3data_to_check(m,18) = 3;
    end
    
    % transpose the data to check
    size3data_to_checkT = size3data_to_check(:,1:24)';
    %toc   
    %mex checkingconstraints_new.cpp       
    size3outcome = checkingconstraints_size3_sorted(size3data_to_checkT,10000);
    size3startindex1 = min(find(size3data_to_check_all(:,1)==0));
    % startindex2 = min(find(outcome_to_save==0));
    size3startindex2 = (size3number_of_generation-1)*10000 + 1;
    size3data_to_check_all(size3startindex1:size3startindex1+10000-1,:) = size3data_to_check;
    size3outcome_to_save(size3startindex2:(size3startindex2+10000-1)) = size3outcome;
    %sum(outcome)
    end
    
    syndatao_size3_old = size3data_to_check_all(find(size3outcome_to_save==1),:);
    extras_size3_old2 = size3data_to_check_all(find(size3outcome_to_save==0),:);
   
    syndatao_size3 = syndatao_size3_old(1:ACSsize3_count,:);
    size3exceeding_index = syndatao_size3(ACSsize3_count,1);

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
    
    %% size4
    size4data_to_check_all = zeros(5000000,8*4+1+4);
    size4outcome_to_save = zeros(5000000,1);
    size4number_of_generation = 0;
    while (sum(size4outcome_to_save)<=ACSsize4_count)
    size4number_of_generation = size4number_of_generation + 1;
    %tic;
    size4data_to_check = zeros(10000,8*4+1+4);
    for m=1:10000
        pi_size4 = pi.*lambda2(:,3);
        pi_size4_renorm = pi_size4./sum(pi_size4);
        hhindexh = randomsample(pi_size4_renorm,rand);
        size4data_to_check(m,33) = hhindexh;
        syn = zeros(4,p+1);
        for hh = 1:4
            memberindexhh = randomsample(w(hhindexh,:),rand);
            size4data_to_check(m,33+hh) = memberindexhh;
            % generating individual level data
            for j = 1:p
                phimj = reshape(phi((hhindexh-1)*L+memberindexhh,j,1:d(j)),1,d(j));
                syn(hh,j) = randomsample(phimj,rand);
            end
        end
        % generating household level data
        lambda1mj = lambda1(hhindexh,:);
        syn(:,p+1) = randomsample(lambda1mj,rand);
        syn_sorted = sortrows(syn,5);
        size4data_to_check(m,3:8) = syn_sorted(1,:);
        size4data_to_check(m,11:16) = syn_sorted(2,:);
        size4data_to_check(m,19:24) = syn_sorted(3,:);
        size4data_to_check(m,27:32) = syn_sorted(4,:);
        size4data_to_check(m,1) = m + (size2number_of_generation+size3number_of_generation+size4number_of_generation)*10000;
        size4data_to_check(m,9) = m + (size2number_of_generation+size3number_of_generation+size4number_of_generation)*10000;
        size4data_to_check(m,17) = m + (size2number_of_generation+size3number_of_generation+size4number_of_generation)*10000;
        size4data_to_check(m,25) = m + (size2number_of_generation+size3number_of_generation+size4number_of_generation)*10000;
        size4data_to_check(m,2) = 1;
        size4data_to_check(m,10) = 2;
        size4data_to_check(m,18) = 3;
        size4data_to_check(m,26) = 4;
    end
    
    % transpose the data to check
    size4data_to_checkT = size4data_to_check(:,1:32)';
    %toc   
    %mex checkingconstraints_new.cpp       
    size4outcome1 = checkingconstraints_size4_sorted_part1(size4data_to_checkT,10000);
    size4outcome2 = checkingconstraints_size4_sorted_part2(size4data_to_checkT,10000);
    size4outcome3 = checkingconstraints_size4_sorted_part3(size4data_to_checkT,10000);
    size4outcome = size4outcome1+size4outcome2+size4outcome3;
    size4startindex1 = min(find(size4data_to_check_all(:,1)==0));
    % startindex2 = min(find(outcome_to_save==0));
    size4startindex2 = (size4number_of_generation-1)*10000 + 1;
    size4data_to_check_all(size4startindex1:size4startindex1+10000-1,:) = size4data_to_check;
    size4outcome_to_save(size4startindex2:(size4startindex2+10000-1)) = size4outcome;
    %sum(outcome)
    end
    
    syndatao_size4_old = size4data_to_check_all(find(size4outcome_to_save==1),:);
    extras_size4_old2 = size4data_to_check_all(find(size4outcome_to_save==0),:);
   
    syndatao_size4 = syndatao_size4_old(1:ACSsize4_count,:);
    size4exceeding_index = syndatao_size4(ACSsize4_count,1);
    
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
    if ((i==4960)||(i==4970)||(i==4980)||(i==4990)||(i==5000)||(i==9960)||(i==9970)||(i==9980)||(i==9990)||(i==10000))
        syndatao_size2_new = zeros(2*ACSsize2_count,8);
        syndatao_size3_new = zeros(3*ACSsize3_count,8);
        syndatao_size4_new = zeros(4*ACSsize4_count,8);
        
        for h=1:ACSsize2_count
            syndatao_size2_new(2*(h-1)+1,:) = syndatao_size2(h,1:8);
            syndatao_size2_new(2*h,:) = syndatao_size2(h,9:16);
        end
        
        for h=1:ACSsize3_count
            syndatao_size3_new(3*(h-1)+1,:) = syndatao_size3(h,1:8);
            syndatao_size3_new(3*(h-1)+2,:) = syndatao_size3(h,9:16);
            syndatao_size3_new(3*h,:) = syndatao_size3(h,17:24);
        end
        for h=1:ACSsize4_count
            syndatao_size4_new(4*(h-1)+1,:) = syndatao_size4(h,1:8);
            syndatao_size4_new(4*(h-1)+2,:) = syndatao_size4(h,9:16);
            syndatao_size4_new(4*(h-1)+3,:) = syndatao_size4(h,17:24);
            syndatao_size4_new(4*h,:) = syndatao_size4(h,25:32);
        end
    end
    
    if i==4960
        syndatao1 = [syndatao_size2_new;syndatao_size3_new;syndatao_size4_new];
    elseif i==4970
        syndatao2 = [syndatao_size2_new;syndatao_size3_new;syndatao_size4_new];
    elseif i==4980
        syndatao3 = [syndatao_size2_new;syndatao_size3_new;syndatao_size4_new];
    elseif i==4990
        syndatao4 = [syndatao_size2_new;syndatao_size3_new;syndatao_size4_new];
    elseif i==5000
        syndatao5 = [syndatao_size2_new;syndatao_size3_new;syndatao_size4_new];
    elseif i==9960
        syndatao6 = [syndatao_size2_new;syndatao_size3_new;syndatao_size4_new];
    elseif i==9970
        syndatao7 = [syndatao_size2_new;syndatao_size3_new;syndatao_size4_new];
    elseif i==9980
        syndatao8 = [syndatao_size2_new;syndatao_size3_new;syndatao_size4_new];
    elseif i==9990
        syndatao9 = [syndatao_size2_new;syndatao_size3_new;syndatao_size4_new];
    elseif i==10000
        syndatao10 = [syndatao_size2_new;syndatao_size3_new;syndatao_size4_new];
    end

    %truesyndatao = syndatao(1:min(find(syndatao(:,1)==0))-1,:);
    %toc
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

HHindex = unique(ACSdatasample(:,1)); 
SSorig = hist(ACSdatasample(:,1),HHindex); 

%profile viewer
%profile off

save('ACSsimulation_size2size3size4_10000HHs_060715output.mat','p','n_sout','nout','piout','K','L','alphaout','syndatao1','syndatao2','syndatao3',...
    'syndatao4','syndatao5','syndatao6','syndatao7','syndatao8','syndatao9','syndatao10',...
    'betaout','wout','data','origdata','SS','piout','z_HH_save','z_member_save','newphiout','HHindex','SSorig','ACSdatasample','elapsed_time','lambda1out','lambda2out',...
    'size2extrasize','size3extrasize','size4extrasize','ACSsize2_count','ACSsize3_count','ACSsize4_count');
toc    
