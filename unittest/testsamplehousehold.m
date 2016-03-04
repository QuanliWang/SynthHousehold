clear all
addpath dist
load last
lambda1 = lambda{1};
lambda2 = lambda{2};

rn = rand(10000,hh_size + 2 + hh_size * p);

tic
    L = size(w,2);
    p = length(d);
    data_to_check = zeros(howmany,8*hh_size+1+hh_size);
    
    %no need to normalize it as it is cared in randomsample_new function
    hhindexh = randomsample_new(pi.*lambda2(:,hh_size-1),rn(:,1));
    data_to_check(:,hh_size * 8 + 1) = hhindexh;
    
    %memberindexhh
    %do random samples for the same probs at the same time
    levelhh = unique(hhindexh); 
    cluster_count = hist(hhindexh,levelhh);
    for i = 1:length(levelhh)
        rntemp = rn(hhindexh == levelhh(i),2:5);
        memberindexhh = randomsample_new(w(levelhh(i),:),rntemp);
        data_to_check(hhindexh == levelhh(i),hh_size * 8 +1 + (1:hh_size))=...
            reshape(memberindexhh,[],hh_size);
        
        % generating household level data
        data_to_check(hhindexh == levelhh(i),(1:hh_size)*8) = ...
            repmat(randomsample_new(lambda1(levelhh(i),:),rn(hhindexh == levelhh(i),6)), 1,hh_size);
        
    end
    
    %phi = reshape(phi,[],p,size(phi,2)); %DIM(phi) = (maxd,p,K*L) or (maxd*p, K*L)
    for hh = 1:hh_size
        group_index = (hhindexh-1)*L+data_to_check(:,hh_size * 8 +1 + hh);
        syn = sampleindividuals(phi,d,group_index,p,rn(:,6 + (hh-1)*p + (1:p)),p); 
        data_to_check(:, (hh-1) * 8 + (3:7)) = syn;
        
        data_to_check(:,(hh-1) * 8 + 1) = (1:howmany) + current_number_of_generation*howmany;
        data_to_check(:,(hh-1) * 8 + 2) = hh;  
        
    end
toc

tic
data_to_check_new = samplehouseholds(phi,w, pi, d, lambda1,lambda2,current_number_of_generation, howmany,hh_size,rn);
toc

temp = data_to_check(:,33:37) - data_to_check_new(:,33:37);
[max(temp(:)) min(temp(:))]

temp = data_to_check(:, (1:hh_size)*8) -data_to_check_new(:, (1:hh_size)*8);
[max(temp(:)) min(temp(:))]

temp = data_to_check(:) -data_to_check_new(:);
[max(temp(:)) min(temp(:))]
