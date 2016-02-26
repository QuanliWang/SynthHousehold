results = load('ACSsimulation_size2size3size4_10000HHs_060715output.mat');

plot(results.alphaout);
plot(results.betaout);

kout = zeros(1,10000);
lout = zeros(1,10000);
for i=1:10000
    kout(i) = length(unique(results.z_HH_save(i,:)));
    lout(i) = length(unique(results.z_member_save(i,:)));
    i
end
plot(kout(5001:10000));
plot(lout(5001:10000));
mean(kout(5001:10000))
quantile(kout(5001:10000),.025)
quantile(kout(5001:10000),.975)

min(lout(5001:10000))
max(lout(5001:10000))

plot(results.n_sout);
mean(results.n_sout(5001:10000))

plot(results.nout);
mean(results.nout(5001:10000))
plot(results.nout(5001:10000));
plot(results.nout(6001:10000));
plot(results.nout(7001:10000));

plot(results.elapsed_time);

origdata = results.origdata(:,1:8);

syndata1 = results.syndatao6;
syndata2 = results.syndatao7;
syndata3 = results.syndatao8;
syndata4 = results.syndatao9;
syndata5 = results.syndatao10;

dlmwrite('ACSsimulation_size234_syndata1_060715.txt',syndata1,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata2_060715.txt',syndata2,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata3_060715.txt',syndata3,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata4_060715.txt',syndata4,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata5_060715.txt',syndata5,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_origdata_060715.txt',origdata,'delimiter',',','precision',5);

dlmwrite('ACSsimulation_size234_nout5001to10000_060715.txt',results.nout(5001:10000),'delimiter',',','precision',5);



%% 8010 - 8100: 10 synthetic datasets
extra_count_8010_8100 = [results.nout_extend(8010),results.nout_extend(8020),results.nout_extend(8030),results.nout_extend(8040),results.nout_extend(8050),...
    results.nout_extend(8060),results.nout_extend(8070),results.nout_extend(8080),results.nout_extend(8090),results.nout_extend(8100)];

dlmwrite('ACSsimulation_size234_syndata1_042214.txt',results.syndatao1,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata2_042214.txt',results.syndatao2,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata3_042214.txt',results.syndatao3,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata4_042214.txt',results.syndatao4,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata5_042214.txt',results.syndatao5,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata6_042214.txt',results.syndatao6,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata7_042214.txt',results.syndatao7,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata8_042214.txt',results.syndatao8,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata9_042214.txt',results.syndatao9,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_syndata10_042214.txt',results.syndatao10,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_size234_origdata_012915.txt',origdata,'delimiter',',','precision',5);

%% 9910 - 10000: 10 synthetic datasets
extra_count_9910_10000 = [results.nout_extend(9910),results.nout_extend(9920),results.nout_extend(9930),results.nout_extend(9940),results.nout_extend(9950),...
    results.nout_extend(9960),results.nout_extend(9970),results.nout_extend(9980),results.nout_extend(9990),results.nout_extend(10000)];

plot(extra_count_8010_8100);
hold on
plot(extra_count_9910_10000);
hold off 
%% check marginals
% sex
origdata_sex = hist(origdata(:,3),[1,2]);
origdata_sex/sum(origdata_sex)

syn9500_sex = hist(syndata9500(:,3),[1,2]);
syn9500_sex/sum(syn9500_sex)

syn9999_sex = hist(syndata9999(:,3),[1,2]);
syn9999_sex/sum(syn9999_sex)

syn10000_sex = hist(syndata10000(:,3),[1,2]);
syn10000_sex/sum(syn10000_sex)

% age
origdata_age = hist(origdata(:,4),1:19);
origdata_age/sum(origdata_age)

syn9500_age = hist(syndata9500(:,4),1:19);
syn9500_age/sum(syn9500_age)

syn9999_age = hist(syndata9999(:,4),1:19);
syn9999_age/sum(syn9999_age)

syn10000_age = hist(syndata10000(:,4),1:19);
syn10000_age/sum(syn10000_age)

% relate
origdata_relate = hist(origdata(:,5),1:12);
origdata_relate/sum(origdata_relate)

syn9500_relate = hist(syndata9500(:,5),1:12);
syn9500_relate/sum(syn9500_relate)

syn9999_relate = hist(syndata9999(:,5),1:12);
syn9999_relate/sum(syn9999_relate)

syn10000_relate = hist(syndata10000(:,5),1:12);
syn10000_relate/sum(syn10000_relate)

% ownership
origdata_own = hist(origdata(:,8),1:2);
origdata_own/sum(origdata_own)

syn9500_own = hist(syndata9500(:,8),1:2);
syn9500_own/sum(syn9500_own)

syn9999_own = hist(syndata9999(:,8),1:2);
syn9999_own/sum(syn9999_own)

syn10000_own = hist(syndata10000(:,8),1:2);
syn10000_own/sum(syn10000_own)

mex size2_recovery.cpp
mex size3_recovery.cpp
mex size4_recovery_part1.cpp
mex size4_recovery_part2.cpp
mex size4_recovery_part3.cpp

%% recovery of combinations
size2_count = 5370;
size3_count = 2504;
size4_count = 2126;

origdata_size2 = origdata(1:size2_count*2,:);
origsize2_outcome = size2_recovery(origdata_size2',size2_count);
origsize2_count = hist(origsize2_outcome,1:11);

syndata1_size2 = syndata1(1:size2_count*2,:);
syndata1size2_outcome = size2_recovery(syndata1_size2',size2_count);
syndata1size2_count = hist(syndata1size2_outcome,1:11);

syndata2_size2 = syndata2(1:size2_count*2,:);
syndata2size2_outcome = size2_recovery(syndata2_size2',size2_count);
syndata2size2_count = hist(syndata2size2_outcome,1:11);

syndata3_size2 = syndata3(1:size2_count*2,:);
syndata3size2_outcome = size2_recovery(syndata3_size2',size2_count);
syndata3size2_count = hist(syndata3size2_outcome,1:11);

syndata4_size2 = syndata4(1:size2_count*2,:);
syndata4size2_outcome = size2_recovery(syndata4_size2',size2_count);
syndata4size2_count = hist(syndata4size2_outcome,1:11);

syndata5_size2 = syndata5(1:size2_count*2,:);
syndata5size2_outcome = size2_recovery(syndata5_size2',size2_count);
syndata5size2_count = hist(syndata5size2_outcome,1:11);

origsize2_count/sum(origsize2_count)

a1=syndata1size2_count/sum(syndata1size2_count)

a2=syndata2size2_count/sum(syndata2size2_count)

a3=syndata3size2_count/sum(syndata3size2_count)

a4=syndata4size2_count/sum(syndata4size2_count)

a5=syndata5size2_count/sum(syndata5size2_count)

(a1(1)+a2(1)+a3(1)+a4(1)+a5(1))/5

(a1(2)+a2(2)+a3(2)+a4(2)+a5(2))/5
%% size 3
origdata_size3 = origdata(size2_count*2+1:size2_count*2+size3_count*3,:);
origsize3_outcome = size3_recovery(origdata_size3',size3_count);
origsize3_count = hist(origsize3_outcome,1:65);

syndata1_size3 = syndata1(size2_count*2+1:size2_count*2+size3_count*3,:);
syndata1size3_outcome = size3_recovery(syndata1_size3',size3_count);
syndata1size3_count = hist(syndata1size3_outcome,1:65);

syndata2_size3 = syndata2(size2_count*2+1:size2_count*2+size3_count*3,:);
syndata2size3_outcome = size3_recovery(syndata2_size3',size3_count);
syndata2size3_count = hist(syndata2size3_outcome,1:65);

syndata3_size3 = syndata3(size2_count*2+1:size2_count*2+size3_count*3,:);
syndata3size3_outcome = size3_recovery(syndata3_size3',size3_count);
syndata3size3_count = hist(syndata3size3_outcome,1:65);

syndata4_size3 = syndata4(size2_count*2+1:size2_count*2+size3_count*3,:);
syndata4size3_outcome = size3_recovery(syndata4_size3',size3_count);
syndata4size3_count = hist(syndata4size3_outcome,1:65);

syndata5_size3 = syndata5(size2_count*2+1:size2_count*2+size3_count*3,:);
syndata5size3_outcome = size3_recovery(syndata5_size3',size3_count);
syndata5size3_count = hist(syndata5size3_outcome,1:65);

origsize3_count/sum(origsize3_count)

b1=syndata1size3_count/sum(syndata1size3_count)

b2=syndata2size3_count/sum(syndata2size3_count)

b3=syndata3size3_count/sum(syndata3size3_count)

b4=syndata4size3_count/sum(syndata4size3_count)

b5=syndata5size3_count/sum(syndata5size3_count)

(b1(1)+b2(1)+b3(1)+b4(1)+b5(1))/5

(b1(11)+b2(11)+b3(11)+b4(11)+b5(11))/5

%% size 4
origdata_size4 = origdata(size2_count*2+size3_count*3+1:size2_count*2+size3_count*3+size4_count*4,:);
origsize4_outcome_part1 = size4_recovery_part1(origdata_size4',size4_count);
origsize4_outcome_part2 = size4_recovery_part2(origdata_size4',size4_count);
origsize4_outcome_part3 = size4_recovery_part3(origdata_size4',size4_count);
origsize4_outcome = origsize4_outcome_part1+origsize4_outcome_part2+origsize4_outcome_part3;
origsize4_count = hist(origsize4_outcome,1:273);

syndata1_size4 = syndata1(size2_count*2+size3_count*3+1:size2_count*2+size3_count*3+size4_count*4,:);
syndata1size4_outcome_part1 = size4_recovery_part1(syndata1_size4',size4_count);
syndata1size4_outcome_part2 = size4_recovery_part2(syndata1_size4',size4_count);
syndata1size4_outcome_part3 = size4_recovery_part3(syndata1_size4',size4_count);
syndata1size4_outcome = syndata1size4_outcome_part1+syndata1size4_outcome_part2+syndata1size4_outcome_part3;
syndata1size4_count = hist(syndata1size4_outcome,1:273);

syndata2_size4 = syndata2(size2_count*2+size3_count*3+1:size2_count*2+size3_count*3+size4_count*4,:);
syndata2size4_outcome_part1 = size4_recovery_part1(syndata2_size4',size4_count);
syndata2size4_outcome_part2 = size4_recovery_part2(syndata2_size4',size4_count);
syndata2size4_outcome_part3 = size4_recovery_part3(syndata2_size4',size4_count);
syndata2size4_outcome = syndata2size4_outcome_part1+syndata2size4_outcome_part2+syndata2size4_outcome_part3;
syndata2size4_count = hist(syndata2size4_outcome,1:273);

syndata3_size4 = syndata3(size2_count*2+size3_count*3+1:size2_count*2+size3_count*3+size4_count*4,:);
syndata3size4_outcome_part1 = size4_recovery_part1(syndata3_size4',size4_count);
syndata3size4_outcome_part2 = size4_recovery_part2(syndata3_size4',size4_count);
syndata3size4_outcome_part3 = size4_recovery_part3(syndata3_size4',size4_count);
syndata3size4_outcome = syndata3size4_outcome_part1+syndata3size4_outcome_part2+syndata3size4_outcome_part3;
syndata3size4_count = hist(syndata3size4_outcome,1:273);

syndata4_size4 = syndata4(size2_count*2+size3_count*3+1:size2_count*2+size3_count*3+size4_count*4,:);
syndata4size4_outcome_part1 = size4_recovery_part1(syndata4_size4',size4_count);
syndata4size4_outcome_part2 = size4_recovery_part2(syndata4_size4',size4_count);
syndata4size4_outcome_part3 = size4_recovery_part3(syndata4_size4',size4_count);
syndata4size4_outcome = syndata4size4_outcome_part1+syndata4size4_outcome_part2+syndata4size4_outcome_part3;
syndata4size4_count = hist(syndata4size4_outcome,1:273);

syndata5_size4 = syndata5(size2_count*2+size3_count*3+1:size2_count*2+size3_count*3+size4_count*4,:);
syndata5size4_outcome_part1 = size4_recovery_part1(syndata5_size4',size4_count);
syndata5size4_outcome_part2 = size4_recovery_part2(syndata5_size4',size4_count);
syndata5size4_outcome_part3 = size4_recovery_part3(syndata5_size4',size4_count);
syndata5size4_outcome = syndata5size4_outcome_part1+syndata5size4_outcome_part2+syndata5size4_outcome_part3;
syndata5size4_count = hist(syndata5size4_outcome,1:273);


origsize4_count/sum(origsize4_count)

c1=syndata1size4_count/sum(syndata1size4_count)

c2=syndata2size4_count/sum(syndata2size4_count)

c3=syndata3size4_count/sum(syndata3size4_count)

c4=syndata4size4_count/sum(syndata4size4_count)

c5=syndata5size4_count/sum(syndata5size4_count)

(c1(1)+c2(1)+c3(1)+c4(1)+c5(1))/5


%%
%SSorig = results.SSorig;
SSorig = results.SS;
SSsyn = [ones(1,5370)*2,ones(1,2504)*3,ones(1,2126)*4];
%SSsyn = [ones(1,5369)*2,ones(1,2503)*3,ones(1,2128)*4];
%SSsyn = SSorig;
%SSsyn = SSorig;
% find one spouse in original and syn1-syn5
one_spouse_orig = zeros(1,10000);
relate_vector_1 = origdata(1:sum(SSorig(1)),7);
if (sum(relate_vector_1==2)==1)
    one_spouse_orig(1) = 1;
end

for h=2:10000
    relate_vector_h = origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7);
    if (sum(relate_vector_h==2)==1)
        one_spouse_orig(h) = 1;
    end
end

one_spouse_syndata1 = zeros(1,10000);
relate_vector_1 = syndata1(1:sum(SSsyn(1)),7);
if (sum(relate_vector_1==2)==1)
    one_spouse_syndata1(1) = 1;
end

for h=2:10000
    relate_vector_h = syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
    if (sum(relate_vector_h==2)==1)
        one_spouse_syndata1(h) = 1;
    end
end

one_spouse_syndata2 = zeros(1,10000);
relate_vector_1 = syndata2(1:sum(SSsyn(1)),7);
if (sum(relate_vector_1==2)==1)
    one_spouse_syndata2(1) = 1;
end

for h=2:10000
    relate_vector_h = syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
    if (sum(relate_vector_h==2)==1)
        one_spouse_syndata2(h) = 1;
    end
end


one_spouse_syndata3 = zeros(1,10000);
relate_vector_1 = syndata3(1:sum(SSsyn(1)),7);
if (sum(relate_vector_1==2)==1)
    one_spouse_syndata3(1) = 1;
end

for h=2:10000
    relate_vector_h = syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
    if (sum(relate_vector_h==2)==1)
        one_spouse_syndata3(h) = 1;
    end
end

one_spouse_syndata4 = zeros(1,10000);
relate_vector_1 = syndata4(1:sum(SSsyn(1)),7);
if (sum(relate_vector_1==2)==1)
    one_spouse_syndata4(1) = 1;
end

for h=2:10000
    relate_vector_h = syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
    if (sum(relate_vector_h==2)==1)
        one_spouse_syndata4(h) = 1;
    end
end


one_spouse_syndata5 = zeros(1,10000);
relate_vector_1 = syndata5(1:sum(SSsyn(1)),7);
if (sum(relate_vector_1==2)==1)
    one_spouse_syndata5(1) = 1;
end

for h=2:10000
    relate_vector_h = syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
    if (sum(relate_vector_h==2)==1)
        one_spouse_syndata5(h) = 1;
    end
end


sum(one_spouse_orig)/10000
[sum(one_spouse_syndata1)/10000,sum(one_spouse_syndata2)/10000,sum(one_spouse_syndata3)/10000,sum(one_spouse_syndata4)/10000,sum(one_spouse_syndata5)/10000]

%% find more than one parents in original and syn1=syn5
morethanone_parent_orig = zeros(1,10000);
relate_vector_1 = origdata(1:sum(SSorig(1)),7);
if ((sum(relate_vector_1==5)+sum(relate_vector_1==6))>=1)
    morethanone_parent_orig(1) = 1;
end

for h=2:10000
    relate_vector_h = origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7);
    if ((sum(relate_vector_h==5)+sum(relate_vector_h==6))>=1)
        morethanone_parent_orig(h) = 1;
    end
end


morethanone_parent_syndata1 = zeros(1,10000);
relate_vector_1 = syndata1(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==5)+sum(relate_vector_1==6))>=1)
    morethanone_parent_syndata1(1) = 1;
end

for h=2:10000
    relate_vector_h = syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
    if ((sum(relate_vector_h==5)+sum(relate_vector_h==6))>=1)
        morethanone_parent_syndata1(h) = 1;
    end
end


morethanone_parent_syndata2 = zeros(1,10000);
relate_vector_1 = syndata2(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==5)+sum(relate_vector_1==6))>=1)
    morethanone_parent_syndata2(1) = 1;
end

for h=2:10000
    relate_vector_h = syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
    if ((sum(relate_vector_h==5)+sum(relate_vector_h==6))>=1)
        morethanone_parent_syndata2(h) = 1;
    end
end



morethanone_parent_syndata3 = zeros(1,10000);
relate_vector_1 = syndata3(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==5)+sum(relate_vector_1==6))>=1)
    morethanone_parent_syndata3(1) = 1;
end

for h=2:10000
    relate_vector_h = syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
    if ((sum(relate_vector_h==5)+sum(relate_vector_h==6))>=1)
        morethanone_parent_syndata3(h) = 1;
    end
end


morethanone_parent_syndata4 = zeros(1,10000);
relate_vector_1 = syndata4(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==5)+sum(relate_vector_1==6))>=1)
    morethanone_parent_syndata4(1) = 1;
end

for h=2:10000
    relate_vector_h = syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
    if ((sum(relate_vector_h==5)+sum(relate_vector_h==6))>=1)
        morethanone_parent_syndata4(h) = 1;
    end
end


morethanone_parent_syndata5 = zeros(1,10000);
relate_vector_1 = syndata5(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==5)+sum(relate_vector_1==6))>=1)
    morethanone_parent_syndata5(1) = 1;
end

for h=2:10000
    relate_vector_h = syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
    if ((sum(relate_vector_h==5)+sum(relate_vector_h==6))>=1)
        morethanone_parent_syndata5(h) = 1;
    end
end

sum(morethanone_parent_orig)/10000
[sum(morethanone_parent_syndata1)/10000,sum(morethanone_parent_syndata2)/10000,sum(morethanone_parent_syndata3)/10000,sum(morethanone_parent_syndata4)/10000,sum(morethanone_parent_syndata5)/10000]


%% find more than one child in original and syn1=syn5
morethanone_child_orig = zeros(1,10000);
relate_vector_1 = origdata(1:sum(SSorig(1)),7);
if ((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)
  morethanone_child_orig(1) = 1;
end

for h=2:10000
relate_vector_h = origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7);
if ((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)
  morethanone_child_orig(h) = 1;
end
end

morethanone_child_syndata1 = zeros(1,10000);
relate_vector_1 = syndata1(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)
  morethanone_child_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)
  morethanone_child_syndata1(h) = 1;
end
end


morethanone_child_syndata2 = zeros(1,10000);
relate_vector_1 = syndata2(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)
  morethanone_child_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)
  morethanone_child_syndata2(h) = 1;
end
end


morethanone_child_syndata3 = zeros(1,10000);
relate_vector_1 = syndata3(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)
  morethanone_child_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)
  morethanone_child_syndata3(h) = 1;
end
end

morethanone_child_syndata4 = zeros(1,10000);
relate_vector_1 = syndata4(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)
  morethanone_child_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)
  morethanone_child_syndata4(h) = 1;
end
end


morethanone_child_syndata5 = zeros(1,10000);
relate_vector_1 = syndata5(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)
  morethanone_child_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)
  morethanone_child_syndata5(h) = 1;
end
end

sum(morethanone_child_orig)/10000
[sum(morethanone_child_syndata1)/10000,sum(morethanone_child_syndata2)/10000,sum(morethanone_child_syndata3)/10000,sum(morethanone_child_syndata4)/10000,sum(morethanone_child_syndata5)/10000]


%% find more than one sbiling in original and syn1=syn5
morethanone_sbiling_orig = zeros(1,10000);
relate_vector_1 = origdata(1:sum(SSorig(1)),7);
if ((sum(relate_vector_1==7)+sum(relate_vector_1==8))>=1)
  morethanone_sbiling_orig(1) = 1;
end

for h=2:10000
relate_vector_h = origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7);
if ((sum(relate_vector_h==7)+sum(relate_vector_h==8))>=1)
  morethanone_sbiling_orig(h) = 1;
end
end

morethanone_sbiling_syndata1 = zeros(1,10000);
relate_vector_1 = syndata1(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==7)+sum(relate_vector_1==8))>=1)
  morethanone_sbiling_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==7)+sum(relate_vector_h==8))>=1)
  morethanone_sbiling_syndata1(h) = 1;
end
end

morethanone_sbiling_syndata2 = zeros(1,10000);
relate_vector_1 = syndata2(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==7)+sum(relate_vector_1==8))>=1)
  morethanone_sbiling_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==7)+sum(relate_vector_h==8))>=1)
  morethanone_sbiling_syndata2(h) = 1;
end
end


morethanone_sbiling_syndata3 = zeros(1,10000);
relate_vector_1 = syndata3(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==7)+sum(relate_vector_1==8))>=1)
  morethanone_sbiling_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==7)+sum(relate_vector_h==8))>=1)
  morethanone_sbiling_syndata3(h) = 1;
end
end


morethanone_sbiling_syndata4 = zeros(1,10000);
relate_vector_1 = syndata4(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==7)+sum(relate_vector_1==8))>=1)
  morethanone_sbiling_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==7)+sum(relate_vector_h==8))>=1)
  morethanone_sbiling_syndata4(h) = 1;
end
end


morethanone_sbiling_syndata5 = zeros(1,10000);
relate_vector_1 = syndata5(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==7)+sum(relate_vector_1==8))>=1)
  morethanone_sbiling_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==7)+sum(relate_vector_h==8))>=1)
  morethanone_sbiling_syndata5(h) = 1;
end
end

sum(morethanone_sbiling_orig)/10000
[sum(morethanone_sbiling_syndata1)/10000,sum(morethanone_sbiling_syndata2)/10000,sum(morethanone_sbiling_syndata3)/10000,sum(morethanone_sbiling_syndata4)/10000,sum(morethanone_sbiling_syndata5)/10000]

%% find morethanone grandchild in original and syn1-syn5
morethanone_grandchild_orig = zeros(1,10000);
relate_vector_1 = origdata(1:sum(SSorig(1)),7);
if (sum(relate_vector_1==9)>=1)
  morethanone_grandchild_orig(1) = 1;
end

for h=2:10000
relate_vector_h = origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7);
if (sum(relate_vector_h==9)>=1)
  morethanone_grandchild_orig(h) = 1;
end
end

sum(morethanone_grandchild_orig);

morethanone_grandchild_syndata1 = zeros(1,10000);
relate_vector_1 = syndata1(1:sum(SSsyn(1)),7);
if (sum(relate_vector_1==9)>=1)
  morethanone_grandchild_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if (sum(relate_vector_h==9)>=1)
  morethanone_grandchild_syndata1(h) = 1;
end
end

morethanone_grandchild_syndata2 = zeros(1,10000);
relate_vector_1 = syndata2(1:sum(SSsyn(1)),7);
if (sum(relate_vector_1==9)>=1)
  morethanone_grandchild_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if (sum(relate_vector_h==9)>=1)
  morethanone_grandchild_syndata2(h) = 1;
end
end


morethanone_grandchild_syndata3 = zeros(1,10000);
relate_vector_1 = syndata3(1:sum(SSsyn(1)),7);
if (sum(relate_vector_1==9)>=1)
  morethanone_grandchild_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if (sum(relate_vector_h==9)>=1)
  morethanone_grandchild_syndata3(h) = 1;
end
end


morethanone_grandchild_syndata4 = zeros(1,10000);
relate_vector_1 = syndata4(1:sum(SSsyn(1)),7);
if (sum(relate_vector_1==9)>=1)
  morethanone_grandchild_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if (sum(relate_vector_h==9)>=1)
  morethanone_grandchild_syndata4(h) = 1;
end
end


morethanone_grandchild_syndata5 = zeros(1,10000);
relate_vector_1 = syndata5(1:sum(SSsyn(1)),7);
if (sum(relate_vector_1==9)>=1)
  morethanone_grandchild_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if (sum(relate_vector_h==9)>=1)
  morethanone_grandchild_syndata5(h) = 1;
end
end

sum(morethanone_grandchild_orig)/10000
[sum(morethanone_grandchild_syndata1)/10000,sum(morethanone_grandchild_syndata2)/10000,sum(morethanone_grandchild_syndata3)/10000,sum(morethanone_grandchild_syndata4)/10000,sum(morethanone_grandchild_syndata5)/10000]


%% find three-generation in original and syn1-syn5
threegen_orig = zeros(1,10000);
relate_vector_1 = origdata(1:sum(SSorig(1)),7);
if ((sum(relate_vector_1==3)+(sum(relate_vector_1==4))>=1&&(sum(relate_vector_1==5)+sum(relate_vector_1==6)>=1))||...
      sum(relate_vector_1==9)>=1)
  threegen_orig(1) = 1;
end

for h=2:10000
relate_vector_h = origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7);
if ((sum(relate_vector_h==3)+(sum(relate_vector_h==4))>=1&&(sum(relate_vector_h==5)+sum(relate_vector_h==6)>=1))||...
      sum(relate_vector_h==9)>=1)
  threegen_orig(h) = 1;
end
end

threegen_syndata1 = zeros(1,10000);
relate_vector_1 = syndata1(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==3)+(sum(relate_vector_1==4))>=1&&(sum(relate_vector_1==5)+sum(relate_vector_1==6)>=1))||...
    sum(relate_vector_1==9)>=1)
  threegen_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==3)+(sum(relate_vector_h==4))>=1&&(sum(relate_vector_h==5)+sum(relate_vector_h==6)>=1))||...
    sum(relate_vector_h==9)>=1)
  threegen_syndata1(h) = 1;
end
end

threegen_syndata2 = zeros(1,10000);
relate_vector_1 = syndata2(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==3)+(sum(relate_vector_1==4))>=1&&(sum(relate_vector_1==5)+sum(relate_vector_1==6)>=1))||...
    sum(relate_vector_1==9)>=1)
  threegen_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==3)+(sum(relate_vector_h==4))>=1&&(sum(relate_vector_h==5)+sum(relate_vector_h==6)>=1))||...
    sum(relate_vector_h==9)>=1)
  threegen_syndata2(h) = 1;
end
end

threegen_syndata3 = zeros(1,10000);
relate_vector_1 = syndata3(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==3)+(sum(relate_vector_1==4))>=1&&(sum(relate_vector_1==5)+sum(relate_vector_1==6)>=1))||...
    sum(relate_vector_1==9)>=1)
  threegen_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==3)+(sum(relate_vector_h==4))>=1&&(sum(relate_vector_h==5)+sum(relate_vector_h==6)>=1))||...
    sum(relate_vector_h==9)>=1)
  threegen_syndata3(h) = 1;
end
end

threegen_syndata4 = zeros(1,10000);
relate_vector_1 = syndata4(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==3)+(sum(relate_vector_1==4))>=1&&(sum(relate_vector_1==5)+sum(relate_vector_1==6)>=1))||...
    sum(relate_vector_1==9)>=1)
  threegen_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==3)+(sum(relate_vector_h==4))>=1&&(sum(relate_vector_h==5)+sum(relate_vector_h==6)>=1))||...
    sum(relate_vector_h==9)>=1)
  threegen_syndata4(h) = 1;
end
end

threegen_syndata5 = zeros(1,10000);
relate_vector_1 = syndata5(1:sum(SSsyn(1)),7);
if ((sum(relate_vector_1==3)+(sum(relate_vector_1==4))>=1&&(sum(relate_vector_1==5)+sum(relate_vector_1==6)>=1))||...
    sum(relate_vector_1==9)>=1)
  threegen_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if ((sum(relate_vector_h==3)+(sum(relate_vector_h==4))>=1&&(sum(relate_vector_h==5)+sum(relate_vector_h==6)>=1))||...
    sum(relate_vector_h==9)>=1)
  threegen_syndata5(h) = 1;
end
end

sum(threegen_orig)/10000
[sum(threegen_syndata1)/10000,sum(threegen_syndata2)/10000,sum(threegen_syndata3)/10000,sum(threegen_syndata4)/10000,sum(threegen_syndata5)/10000]


%% find one spouseHHHwhite in original and syn1-syn5
one_spouseHHHwhite_orig = zeros(1,10000);
relate_vector_1 = [origdata(1:sum(SSorig(1)),7),origdata(1:sum(SSorig(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==1)
  one_spouseHHHwhite_orig(1) = 1;
end

for h=2:10000
relate_vector_h = [origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==1)
  one_spouseHHHwhite_orig(h) = 1;
end
end

one_spouseHHHwhite_syndata1 = zeros(1,10000);
relate_vector_1 = [syndata1(1:sum(SSsyn(1)),7),syndata1(1:sum(SSsyn(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==1)
  one_spouseHHHwhite_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==1)
  one_spouseHHHwhite_syndata1(h) = 1;
end
end

one_spouseHHHwhite_syndata2 = zeros(1,10000);
relate_vector_1 = [syndata2(1:sum(SSsyn(1)),7),syndata2(1:sum(SSsyn(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==1)
  one_spouseHHHwhite_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==1)
  one_spouseHHHwhite_syndata2(h) = 1;
end
end

one_spouseHHHwhite_syndata3 = zeros(1,10000);
relate_vector_1 = [syndata3(1:sum(SSsyn(1)),7),syndata3(1:sum(SSsyn(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==1)
  one_spouseHHHwhite_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==1)
  one_spouseHHHwhite_syndata3(h) = 1;
end
end

one_spouseHHHwhite_syndata4 = zeros(1,10000);
relate_vector_1 = [syndata4(1:sum(SSsyn(1)),7),syndata4(1:sum(SSsyn(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==1)
  one_spouseHHHwhite_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==1)
  one_spouseHHHwhite_syndata4(h) = 1;
end
end

one_spouseHHHwhite_syndata5 = zeros(1,10000);
relate_vector_1 = [syndata5(1:sum(SSsyn(1)),7),syndata5(1:sum(SSsyn(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==1)
  one_spouseHHHwhite_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==1)
  one_spouseHHHwhite_syndata5(h) = 1;
end
end

sum(one_spouseHHHwhite_orig)/10000
[sum(one_spouseHHHwhite_syndata1)/10000,sum(one_spouseHHHwhite_syndata2)/10000,sum(one_spouseHHHwhite_syndata3)/10000,sum(one_spouseHHHwhite_syndata4)/10000,sum(one_spouseHHHwhite_syndata5)/10000]

%% find one spouseHHHblack in original and syn1-syn5
one_spouseHHHblack_orig = zeros(1,10000);
relate_vector_1 = [origdata(1:sum(SSorig(1)),7),origdata(1:sum(SSorig(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==2)
  one_spouseHHHblack_orig(1) = 1;
end

for h=2:10000
relate_vector_h = [origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==2)
  one_spouseHHHblack_orig(h) = 1;
end
end

one_spouseHHHblack_syndata1 = zeros(1,10000);
relate_vector_1 = [syndata1(1:sum(SSsyn(1)),7),syndata1(1:sum(SSsyn(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==2)
  one_spouseHHHblack_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==2)
  one_spouseHHHblack_syndata1(h) = 1;
end
end

one_spouseHHHblack_syndata2 = zeros(1,10000);
relate_vector_1 = [syndata2(1:sum(SSsyn(1)),7),syndata2(1:sum(SSsyn(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==2)
  one_spouseHHHblack_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==2)
  one_spouseHHHblack_syndata2(h) = 1;
end
end

one_spouseHHHblack_syndata3 = zeros(1,10000);
relate_vector_1 = [syndata3(1:sum(SSsyn(1)),7),syndata3(1:sum(SSsyn(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==2)
  one_spouseHHHblack_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==2)
  one_spouseHHHblack_syndata3(h) = 1;
end
end

one_spouseHHHblack_syndata4 = zeros(1,10000);
relate_vector_1 = [syndata4(1:sum(SSsyn(1)),7),syndata4(1:sum(SSsyn(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==2)
  one_spouseHHHblack_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==2)
  one_spouseHHHblack_syndata4(h) = 1;
end
end

one_spouseHHHblack_syndata5 = zeros(1,10000);
relate_vector_1 = [syndata5(1:sum(SSsyn(1)),7),syndata5(1:sum(SSsyn(1)),4)];
if (sum(relate_vector_1(:,1)==2)==1&&relate_vector_1(1,2)==2)
  one_spouseHHHblack_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (sum(relate_vector_h(:,1)==2)==1&&relate_vector_h(1,2)==2)
  one_spouseHHHblack_syndata5(h) = 1;
end
end

sum(one_spouseHHHblack_orig)/10000
[sum(one_spouseHHHblack_syndata1)/10000,sum(one_spouseHHHblack_syndata2)/10000,sum(one_spouseHHHblack_syndata3)/10000,sum(one_spouseHHHblack_syndata4)/10000,sum(one_spouseHHHblack_syndata5)/10000]

%% find one bothcouplewhite in original and syn1-syn5
one_bothcouplewhite_orig = zeros(1,10000);
relate_vector_1 = [origdata(1:sum(SSorig(1)),7),origdata(1:sum(SSorig(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1)
  one_bothcouplewhite_orig(1) = 1;
end

for h=2:10000
relate_vector_h = [origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1)
  one_bothcouplewhite_orig(h) = 1;
end
end

one_bothcouplewhite_syndata1 = zeros(1,10000);
relate_vector_1 = [syndata1(1:sum(SSsyn(1)),7),syndata1(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1)
  one_bothcouplewhite_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1)
  one_bothcouplewhite_syndata1(h) = 1;
end
end

one_bothcouplewhite_syndata2 = zeros(1,10000);
relate_vector_1 = [syndata2(1:sum(SSsyn(1)),7),syndata2(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1)
  one_bothcouplewhite_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1)
  one_bothcouplewhite_syndata2(h) = 1;
end
end

one_bothcouplewhite_syndata3 = zeros(1,10000);
relate_vector_1 = [syndata3(1:sum(SSsyn(1)),7),syndata3(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1)
  one_bothcouplewhite_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1)
  one_bothcouplewhite_syndata3(h) = 1;
end
end

one_bothcouplewhite_syndata4 = zeros(1,10000);
relate_vector_1 = [syndata4(1:sum(SSsyn(1)),7),syndata4(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1)
  one_bothcouplewhite_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1)
  one_bothcouplewhite_syndata4(h) = 1;
end
end

one_bothcouplewhite_syndata5 = zeros(1,10000);
relate_vector_1 = [syndata5(1:sum(SSsyn(1)),7),syndata5(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1)
  one_bothcouplewhite_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1)
  one_bothcouplewhite_syndata5(h) = 1;
end
end

sum(one_bothcouplewhite_orig)/10000
[sum(one_bothcouplewhite_syndata1)/10000,sum(one_bothcouplewhite_syndata2)/10000,sum(one_bothcouplewhite_syndata3)/10000,sum(one_bothcouplewhite_syndata4)/10000,sum(one_bothcouplewhite_syndata5)/10000]


%% find one sameracecouple in original and syn1-syn5
one_sameracecouple_orig = zeros(1,10000);
relate_vector_1 = [origdata(1:sum(SSorig(1)),7),origdata(1:sum(SSorig(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==relate_vector_1(2,2))
  one_sameracecouple_orig(1) = 1;
end

for h=2:10000
relate_vector_h = [origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==relate_vector_h(2,2))
  one_sameracecouple_orig(h) = 1;
end
end

one_sameracecouple_syndata1 = zeros(1,10000);
relate_vector_1 = [syndata1(1:sum(SSsyn(1)),7),syndata1(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==relate_vector_1(2,2))
  one_sameracecouple_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==relate_vector_h(2,2))
  one_sameracecouple_syndata1(h) = 1;
end
end

one_sameracecouple_syndata2 = zeros(1,10000);
relate_vector_1 = [syndata2(1:sum(SSsyn(1)),7),syndata2(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==relate_vector_1(2,2))
  one_sameracecouple_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==relate_vector_h(2,2))
  one_sameracecouple_syndata2(h) = 1;
end
end

one_sameracecouple_syndata3 = zeros(1,10000);
relate_vector_1 = [syndata3(1:sum(SSsyn(1)),7),syndata3(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==relate_vector_1(2,2))
  one_sameracecouple_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==relate_vector_h(2,2))
  one_sameracecouple_syndata3(h) = 1;
end
end

one_sameracecouple_syndata4 = zeros(1,10000);
relate_vector_1 = [syndata4(1:sum(SSsyn(1)),7),syndata4(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==relate_vector_1(2,2))
  one_sameracecouple_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==relate_vector_h(2,2))
  one_sameracecouple_syndata4(h) = 1;
end
end

one_sameracecouple_syndata5 = zeros(1,10000);
relate_vector_1 = [syndata5(1:sum(SSsyn(1)),7),syndata5(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==relate_vector_1(2,2))
  one_sameracecouple_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==relate_vector_h(2,2))
  one_sameracecouple_syndata5(h) = 1;
end
end

sum(one_sameracecouple_orig)/10000
[sum(one_sameracecouple_syndata1)/10000,sum(one_sameracecouple_syndata2)/10000,sum(one_sameracecouple_syndata3)/10000,sum(one_sameracecouple_syndata4)/10000,sum(one_sameracecouple_syndata5)/10000]

%% find one whiteandnonewhitecouple in original and syn1-syn5
one_whiteandnonewhitecouple_orig = zeros(1,10000);
relate_vector_1 = [origdata(1:sum(SSorig(1)),7),origdata(1:sum(SSorig(1)),4)];
if (relate_vector_1(2,1)==2&&((relate_vector_1(1,2)==1&&relate_vector_1(2,2)~=1)||(relate_vector_1(1,2)~=1&&relate_vector_1(2,2)==1)))
  one_whiteandnonewhitecouple_orig(1) = 1;
end

for h=2:10000
relate_vector_h = [origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),4)];
if (relate_vector_h(2,1)==2&&((relate_vector_h(1,2)==1&&relate_vector_h(2,2)~=1)||(relate_vector_h(1,2)~=1&&relate_vector_h(2,2)==1)))
  one_whiteandnonewhitecouple_orig(h) = 1;
end
end

one_whiteandnonewhitecouple_syndata1 = zeros(1,10000);
relate_vector_1 = [syndata1(1:sum(SSsyn(1)),7),syndata1(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&((relate_vector_1(1,2)==1&&relate_vector_1(2,2)~=1)||(relate_vector_1(1,2)~=1&&relate_vector_1(2,2)==1)))
  one_whiteandnonewhitecouple_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&((relate_vector_h(1,2)==1&&relate_vector_h(2,2)~=1)||(relate_vector_h(1,2)~=1&&relate_vector_h(2,2)==1)))
  one_whiteandnonewhitecouple_syndata1(h) = 1;
end
end

one_whiteandnonewhitecouple_syndata2 = zeros(1,10000);
relate_vector_1 = [syndata2(1:sum(SSsyn(1)),7),syndata2(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&((relate_vector_1(1,2)==1&&relate_vector_1(2,2)~=1)||(relate_vector_1(1,2)~=1&&relate_vector_1(2,2)==1)))
  one_whiteandnonewhitecouple_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&((relate_vector_h(1,2)==1&&relate_vector_h(2,2)~=1)||(relate_vector_h(1,2)~=1&&relate_vector_h(2,2)==1)))
  one_whiteandnonewhitecouple_syndata2(h) = 1;
end
end

one_whiteandnonewhitecouple_syndata3 = zeros(1,10000);
relate_vector_1 = [syndata3(1:sum(SSsyn(1)),7),syndata3(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&((relate_vector_1(1,2)==1&&relate_vector_1(2,2)~=1)||(relate_vector_1(1,2)~=1&&relate_vector_1(2,2)==1)))
  one_whiteandnonewhitecouple_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&((relate_vector_h(1,2)==1&&relate_vector_h(2,2)~=1)||(relate_vector_h(1,2)~=1&&relate_vector_h(2,2)==1)))
  one_whiteandnonewhitecouple_syndata3(h) = 1;
end
end

one_whiteandnonewhitecouple_syndata4 = zeros(1,10000);
relate_vector_1 = [syndata4(1:sum(SSsyn(1)),7),syndata4(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&((relate_vector_1(1,2)==1&&relate_vector_1(2,2)~=1)||(relate_vector_1(1,2)~=1&&relate_vector_1(2,2)==1)))
  one_whiteandnonewhitecouple_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&((relate_vector_h(1,2)==1&&relate_vector_h(2,2)~=1)||(relate_vector_h(1,2)~=1&&relate_vector_h(2,2)==1)))
  one_whiteandnonewhitecouple_syndata4(h) = 1;
end
end

one_whiteandnonewhitecouple_syndata5 = zeros(1,10000);
relate_vector_1 = [syndata5(1:sum(SSsyn(1)),7),syndata5(1:sum(SSsyn(1)),4)];
if (relate_vector_1(2,1)==2&&((relate_vector_1(1,2)==1&&relate_vector_1(2,2)~=1)||(relate_vector_1(1,2)~=1&&relate_vector_1(2,2)==1)))
  one_whiteandnonewhitecouple_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4)];
if (relate_vector_h(2,1)==2&&((relate_vector_h(1,2)==1&&relate_vector_h(2,2)~=1)||(relate_vector_h(1,2)~=1&&relate_vector_h(2,2)==1)))
  one_whiteandnonewhitecouple_syndata5(h) = 1;
end
end

sum(one_whiteandnonewhitecouple_orig)/10000
[sum(one_whiteandnonewhitecouple_syndata1)/10000,sum(one_whiteandnonewhitecouple_syndata2)/10000,sum(one_whiteandnonewhitecouple_syndata3)/10000,sum(one_whiteandnonewhitecouple_syndata4)/10000,sum(one_whiteandnonewhitecouple_syndata5)/10000]


%% find one bothcouplewhiteown in original and syn1-syn5
one_bothcouplewhiteown_orig = zeros(1,10000);
relate_vector_1 = [origdata(1:sum(SSorig(1)),7),origdata(1:sum(SSorig(1)),4),origdata(1:sum(SSorig(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1&&relate_vector_1(1,3)==1)
  one_bothcouplewhiteown_orig(1) = 1;
end

for h=2:10000
relate_vector_h = [origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),4),origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1&&relate_vector_h(1,3)==1)
  one_bothcouplewhiteown_orig(h) = 1;
end
end

one_bothcouplewhiteown_syndata1 = zeros(1,10000);
relate_vector_1 = [syndata1(1:sum(SSsyn(1)),7),syndata1(1:sum(SSsyn(1)),4),syndata1(1:sum(SSsyn(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1&&relate_vector_1(1,3)==1)
  one_bothcouplewhiteown_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4),syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1&&relate_vector_h(1,3)==1)
  one_bothcouplewhiteown_syndata1(h) = 1;
end
end

one_bothcouplewhiteown_syndata2 = zeros(1,10000);
relate_vector_1 = [syndata2(1:sum(SSsyn(1)),7),syndata2(1:sum(SSsyn(1)),4),syndata2(1:sum(SSsyn(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1&&relate_vector_1(1,3)==1)
  one_bothcouplewhiteown_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4),syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1&&relate_vector_h(1,3)==1)
  one_bothcouplewhiteown_syndata2(h) = 1;
end
end

one_bothcouplewhiteown_syndata3 = zeros(1,10000);
relate_vector_1 = [syndata3(1:sum(SSsyn(1)),7),syndata3(1:sum(SSsyn(1)),4),syndata3(1:sum(SSsyn(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1&&relate_vector_1(1,3)==1)
  one_bothcouplewhiteown_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4),syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1&&relate_vector_h(1,3)==1)
  one_bothcouplewhiteown_syndata3(h) = 1;
end
end

one_bothcouplewhiteown_syndata4 = zeros(1,10000);
relate_vector_1 = [syndata4(1:sum(SSsyn(1)),7),syndata4(1:sum(SSsyn(1)),4),syndata4(1:sum(SSsyn(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1&&relate_vector_1(1,3)==1)
  one_bothcouplewhiteown_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4),syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1&&relate_vector_h(1,3)==1)
  one_bothcouplewhiteown_syndata4(h) = 1;
end
end

one_bothcouplewhiteown_syndata5 = zeros(1,10000);
relate_vector_1 = [syndata5(1:sum(SSsyn(1)),7),syndata5(1:sum(SSsyn(1)),4),syndata5(1:sum(SSsyn(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)==1&&relate_vector_1(2,2)==1&&relate_vector_1(1,3)==1)
  one_bothcouplewhiteown_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4),syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)==1&&relate_vector_h(2,2)==1&&relate_vector_h(1,3)==1)
  one_bothcouplewhiteown_syndata5(h) = 1;
end
end

sum(one_bothcouplewhiteown_orig)/10000
[sum(one_bothcouplewhiteown_syndata1)/10000,sum(one_bothcouplewhiteown_syndata2)/10000,sum(one_bothcouplewhiteown_syndata3)/10000,sum(one_bothcouplewhiteown_syndata4)/10000,sum(one_bothcouplewhiteown_syndata5)/10000]

%% find one nonwhitecoupleown in original and syn1-syn5
one_nonwhitecoupleown_orig = zeros(1,10000);
relate_vector_1 = [origdata(1:sum(SSorig(1)),7),origdata(1:sum(SSorig(1)),4),origdata(1:sum(SSorig(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)~=1&&relate_vector_1(2,2)~=1&&relate_vector_1(1,3)==1)
  one_nonwhitecoupleown_orig(1) = 1;
end

for h=2:10000
relate_vector_h = [origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),4),origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)~=1&&relate_vector_h(2,2)~=1&&relate_vector_h(1,3)==1)
  one_nonwhitecoupleown_orig(h) = 1;
end
end

one_nonwhitecoupleown_syndata1 = zeros(1,10000);
relate_vector_1 = [syndata1(1:sum(SSsyn(1)),7),syndata1(1:sum(SSsyn(1)),4),syndata1(1:sum(SSsyn(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)~=1&&relate_vector_1(2,2)~=1&&relate_vector_1(1,3)==1)
  one_nonwhitecoupleown_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4),syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)~=1&&relate_vector_h(2,2)~=1&&relate_vector_h(1,3)==1)
  one_nonwhitecoupleown_syndata1(h) = 1;
end
end

one_nonwhitecoupleown_syndata2 = zeros(1,10000);
relate_vector_1 = [syndata2(1:sum(SSsyn(1)),7),syndata2(1:sum(SSsyn(1)),4),syndata2(1:sum(SSsyn(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)~=1&&relate_vector_1(2,2)~=1&&relate_vector_1(1,3)==1)
  one_nonwhitecoupleown_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4),syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)~=1&&relate_vector_h(2,2)~=1&&relate_vector_h(1,3)==1)
  one_nonwhitecoupleown_syndata2(h) = 1;
end
end

one_nonwhitecoupleown_syndata3 = zeros(1,10000);
relate_vector_1 = [syndata3(1:sum(SSsyn(1)),7),syndata3(1:sum(SSsyn(1)),4),syndata3(1:sum(SSsyn(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)~=1&&relate_vector_1(2,2)~=1&&relate_vector_1(1,3)==1)
  one_nonwhitecoupleown_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4),syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)~=1&&relate_vector_h(2,2)~=1&&relate_vector_h(1,3)==1)
  one_nonwhitecoupleown_syndata3(h) = 1;
end
end

one_nonwhitecoupleown_syndata4 = zeros(1,10000);
relate_vector_1 = [syndata4(1:sum(SSsyn(1)),7),syndata4(1:sum(SSsyn(1)),4),syndata4(1:sum(SSsyn(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)~=1&&relate_vector_1(2,2)~=1&&relate_vector_1(1,3)==1)
  one_nonwhitecoupleown_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4),syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)~=1&&relate_vector_h(2,2)~=1&&relate_vector_h(1,3)==1)
  one_nonwhitecoupleown_syndata4(h) = 1;
end
end

one_nonwhitecoupleown_syndata5 = zeros(1,10000);
relate_vector_1 = [syndata5(1:sum(SSsyn(1)),7),syndata5(1:sum(SSsyn(1)),4),syndata5(1:sum(SSsyn(1)),8)];
if (relate_vector_1(2,1)==2&&relate_vector_1(1,2)~=1&&relate_vector_1(2,2)~=1&&relate_vector_1(1,3)==1)
  one_nonwhitecoupleown_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),4),syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),8)];
if (relate_vector_h(2,1)==2&&relate_vector_h(1,2)~=1&&relate_vector_h(2,2)~=1&&relate_vector_h(1,3)==1)
  one_nonwhitecoupleown_syndata5(h) = 1;
end
end

sum(one_nonwhitecoupleown_orig)/10000
[sum(one_nonwhitecoupleown_syndata1)/10000,sum(one_nonwhitecoupleown_syndata2)/10000,sum(one_nonwhitecoupleown_syndata3)/10000,sum(one_nonwhitecoupleown_syndata4)/10000,sum(one_nonwhitecoupleown_syndata5)/10000]

%% find more than one child but no spouse in original and syn1=syn5
nospouse_morethanonechild_orig = zeros(1,10000);
relate_vector_1 = origdata(1:sum(SSorig(1)),7);
if (((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)&&(sum(relate_vector_1==2)==0))
  nospouse_morethanonechild_orig(1) = 1;
end

for h=2:10000
relate_vector_h = origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7);
if (((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)&&(sum(relate_vector_h==2)==0))
  nospouse_morethanonechild_orig(h) = 1;
end
end

nospouse_morethanonechild_syndata1 = zeros(1,10000);
relate_vector_1 = syndata1(1:sum(SSsyn(1)),7);
if (((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)&&(sum(relate_vector_1==2)==0))
  nospouse_morethanonechild_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if (((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)&&(sum(relate_vector_h==2)==0))
  nospouse_morethanonechild_syndata1(h) = 1;
end
end


nospouse_morethanonechild_syndata2 = zeros(1,10000);
relate_vector_1 = syndata2(1:sum(SSsyn(1)),7);
if (((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)&&(sum(relate_vector_1==2)==0))
  nospouse_morethanonechild_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if (((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)&&(sum(relate_vector_h==2)==0))
  nospouse_morethanonechild_syndata2(h) = 1;
end
end



nospouse_morethanonechild_syndata3 = zeros(1,10000);
relate_vector_1 = syndata3(1:sum(SSsyn(1)),7);
if (((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)&&(sum(relate_vector_1==2)==0))
  nospouse_morethanonechild_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if (((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)&&(sum(relate_vector_h==2)==0))
  nospouse_morethanonechild_syndata3(h) = 1;
end
end

nospouse_morethanonechild_syndata4 = zeros(1,10000);
relate_vector_1 = syndata4(1:sum(SSsyn(1)),7);
if (((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)&&(sum(relate_vector_1==2)==0))
  nospouse_morethanonechild_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if (((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)&&(sum(relate_vector_h==2)==0))
  nospouse_morethanonechild_syndata4(h) = 1;
end
end


nospouse_morethanonechild_syndata5 = zeros(1,10000);
relate_vector_1 = syndata5(1:sum(SSsyn(1)),7);
if (((sum(relate_vector_1==3)+sum(relate_vector_1==4))>=1)&&(sum(relate_vector_1==2)==0))
  nospouse_morethanonechild_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7);
if (((sum(relate_vector_h==3)+sum(relate_vector_h==4))>=1)&&(sum(relate_vector_h==2)==0))
  nospouse_morethanonechild_syndata5(h) = 1;
end
end

sum(nospouse_morethanonechild_orig)/10000
[sum(nospouse_morethanonechild_syndata1)/10000,sum(nospouse_morethanonechild_syndata2)/10000,sum(nospouse_morethanonechild_syndata3)/10000,sum(nospouse_morethanonechild_syndata4)/10000,sum(nospouse_morethanonechild_syndata5)/10000]

%% find more than one child but no spouse in original and syn1=syn5
onlymom_morethanonechild_orig = zeros(1,10000);
relate_vector_1 = [origdata(1:sum(SSorig(1)),7),origdata(1:sum(SSorig(1)),3)];
if (((sum(relate_vector_1(:,1)==3)+sum(relate_vector_1(:,1)==4))>=1)&&(sum(relate_vector_1(:,1)==2)==0)&&(relate_vector_1(1,2)==2))
  onlymom_morethanonechild_orig(1) = 1;
end

for h=2:10000
relate_vector_h = [origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),3)];
if (((sum(relate_vector_h(:,1)==3)+sum(relate_vector_h(:,1)==4))>=1)&&(sum(relate_vector_h(:,1)==2)==0)&&(relate_vector_h(1,2)==2))
  onlymom_morethanonechild_orig(h) = 1;
end
end

onlymom_morethanonechild_syndata1 = zeros(1,10000);
relate_vector_1 = [syndata1(1:sum(SSsyn(1)),7),syndata1(1:sum(SSsyn(1)),3)];
if (((sum(relate_vector_1(:,1)==3)+sum(relate_vector_1(:,1)==4))>=1)&&(sum(relate_vector_1(:,1)==2)==0)&&(relate_vector_1(1,2)==2))
  onlymom_morethanonechild_syndata1(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata1(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),3)];
if (((sum(relate_vector_h(:,1)==3)+sum(relate_vector_h(:,1)==4))>=1)&&(sum(relate_vector_h(:,1)==2)==0)&&(relate_vector_h(1,2)==2))
  onlymom_morethanonechild_syndata1(h) = 1;
end
end

onlymom_morethanonechild_syndata2 = zeros(1,10000);
relate_vector_1 = [syndata2(1:sum(SSsyn(1)),7),syndata2(1:sum(SSsyn(1)),3)];
if (((sum(relate_vector_1(:,1)==3)+sum(relate_vector_1(:,1)==4))>=1)&&(sum(relate_vector_1(:,1)==2)==0)&&(relate_vector_1(1,2)==2))
  onlymom_morethanonechild_syndata2(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata2(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),3)];
if (((sum(relate_vector_h(:,1)==3)+sum(relate_vector_h(:,1)==4))>=1)&&(sum(relate_vector_h(:,1)==2)==0)&&(relate_vector_h(1,2)==2))
  onlymom_morethanonechild_syndata2(h) = 1;
end
end

onlymom_morethanonechild_syndata3 = zeros(1,10000);
relate_vector_1 = [syndata3(1:sum(SSsyn(1)),7),syndata3(1:sum(SSsyn(1)),3)];
if (((sum(relate_vector_1(:,1)==3)+sum(relate_vector_1(:,1)==4))>=1)&&(sum(relate_vector_1(:,1)==2)==0)&&(relate_vector_1(1,2)==2))
  onlymom_morethanonechild_syndata3(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata3(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),3)];
if (((sum(relate_vector_h(:,1)==3)+sum(relate_vector_h(:,1)==4))>=1)&&(sum(relate_vector_h(:,1)==2)==0)&&(relate_vector_h(1,2)==2))
  onlymom_morethanonechild_syndata3(h) = 1;
end
end


onlymom_morethanonechild_syndata4 = zeros(1,10000);
relate_vector_1 = [syndata4(1:sum(SSsyn(1)),7),syndata4(1:sum(SSsyn(1)),3)];
if (((sum(relate_vector_1(:,1)==3)+sum(relate_vector_1(:,1)==4))>=1)&&(sum(relate_vector_1(:,1)==2)==0)&&(relate_vector_1(1,2)==2))
  onlymom_morethanonechild_syndata4(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata4(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),3)];
if (((sum(relate_vector_h(:,1)==3)+sum(relate_vector_h(:,1)==4))>=1)&&(sum(relate_vector_h(:,1)==2)==0)&&(relate_vector_h(1,2)==2))
  onlymom_morethanonechild_syndata4(h) = 1;
end
end


onlymom_morethanonechild_syndata5 = zeros(1,10000);
relate_vector_1 = [syndata5(1:sum(SSsyn(1)),7),syndata5(1:sum(SSsyn(1)),3)];
if (((sum(relate_vector_1(:,1)==3)+sum(relate_vector_1(:,1)==4))>=1)&&(sum(relate_vector_1(:,1)==2)==0)&&(relate_vector_1(1,2)==2))
  onlymom_morethanonechild_syndata5(1) = 1;
end

for h=2:10000
relate_vector_h = [syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),7),syndata5(sum(SSsyn(1:h-1))+1:sum(SSsyn(1:h)),3)];
if (((sum(relate_vector_h(:,1)==3)+sum(relate_vector_h(:,1)==4))>=1)&&(sum(relate_vector_h(:,1)==2)==0)&&(relate_vector_h(1,2)==2))
  onlymom_morethanonechild_syndata5(h) = 1;
end
end

sum(onlymom_morethanonechild_orig)/10000
[sum(onlymom_morethanonechild_syndata1)/10000,sum(onlymom_morethanonechild_syndata2)/10000,sum(onlymom_morethanonechild_syndata3)/10000,sum(onlymom_morethanonechild_syndata4)/10000,sum(onlymom_morethanonechild_syndata5)/10000]

%% find one spouseagediff in original and syn1-syn5
%one_spouseagediff_orig_count = zeros(1,10000);
one_spouseagediff_orig = 0;
relate_vector_1 = [origdata(1:sum(SSorig(1)),7),origdata(1:sum(SSorig(1)),6)];
if (relate_vector_1(2,1)==2)
  %one_spouseagediff_orig_count(1) = 1;
  one_spouseagediff_orig = [one_spouseagediff_orig,(relate_vector_1(1,2)-relate_vector_1(2,2))];
end

for h=2:10000
relate_vector_h = [origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),origdata(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),6)];
if (relate_vector_h(2,1)==2)
  one_spouseagediff_orig = [one_spouseagediff_orig,(relate_vector_h(1,2)-relate_vector_h(2,2))];
end
end

%one_spouseagediff_syndata1 = zeros(1,10000);
one_spouseagediff_syndata1 = 0;
relate_vector_1 = [syndata1(1:sum(SSorig(1)),7),syndata1(1:sum(SSorig(1)),6)];
if (relate_vector_1(2,1)==2)
  one_spouseagediff_syndata1 = [one_spouseagediff_syndata1,(relate_vector_1(1,2)-relate_vector_1(2,2))];
end

for h=2:10000
relate_vector_h = [syndata1(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),syndata1(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),6)];
if (relate_vector_h(2,1)==2)
  one_spouseagediff_syndata1 = [one_spouseagediff_syndata1,(relate_vector_h(1,2)-relate_vector_h(2,2))];
end
end

%one_spouseagediff_syndata2 = zeros(1,10000);
one_spouseagediff_syndata2 = 0;
relate_vector_1 = [syndata2(1:sum(SSorig(1)),7),syndata2(1:sum(SSorig(1)),6)];
if (relate_vector_1(2,1)==2)
  one_spouseagediff_syndata2 = [one_spouseagediff_syndata2,(relate_vector_1(1,2)-relate_vector_1(2,2))];
end

for h=2:10000
relate_vector_h = [syndata2(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),syndata2(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),6)];
if (relate_vector_h(2,1)==2)
  one_spouseagediff_syndata2 = [one_spouseagediff_syndata2,(relate_vector_h(1,2)-relate_vector_h(2,2))];
end
end

%one_spouseagediff_syndata3 = zeros(1,10000);
one_spouseagediff_syndata3 = 0;
relate_vector_1 = [syndata3(1:sum(SSorig(1)),7),syndata3(1:sum(SSorig(1)),6)];
if (relate_vector_1(2,1)==2)
  one_spouseagediff_syndata3 = [one_spouseagediff_syndata3,(relate_vector_1(1,2)-relate_vector_1(2,2))];
end

for h=2:10000
relate_vector_h = [syndata3(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),syndata3(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),6)];
if (relate_vector_h(2,1)==2)
  one_spouseagediff_syndata3 = [one_spouseagediff_syndata3,(relate_vector_h(1,2)-relate_vector_h(2,2))];
end
end

%one_spouseagediff_syndata4 = zeros(1,10000);
one_spouseagediff_syndata4 = 0;
relate_vector_1 = [syndata4(1:sum(SSorig(1)),7),syndata4(1:sum(SSorig(1)),6)];
if (relate_vector_1(2,1)==2)
  one_spouseagediff_syndata4 = [one_spouseagediff_syndata4,(relate_vector_1(1,2)-relate_vector_1(2,2))];
end

for h=2:10000
relate_vector_h = [syndata4(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),syndata4(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),6)];
if (relate_vector_h(2,1)==2)
  one_spouseagediff_syndata4 = [one_spouseagediff_syndata4,(relate_vector_h(1,2)-relate_vector_h(2,2))];
end
end

%one_spouseagediff_syndata5 = zeros(1,10000);
one_spouseagediff_syndata5 = 0;
relate_vector_1 = [syndata5(1:sum(SSorig(1)),7),syndata5(1:sum(SSorig(1)),6)];
if (relate_vector_1(2,1)==2)
  one_spouseagediff_syndata5 = [one_spouseagediff_syndata5,(relate_vector_1(1,2)-relate_vector_1(2,2))];
end

for h=2:10000
relate_vector_h = [syndata5(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),7),syndata5(sum(SSorig(1:h-1))+1:sum(SSorig(1:h)),6)];
if (relate_vector_h(2,1)==2)
  one_spouseagediff_syndata5 = [one_spouseagediff_syndata5,(relate_vector_h(1,2)-relate_vector_h(2,2))];
end
end


dlmwrite('DA_spouseagediff_orig.txt',one_spouseagediff_orig','delimiter',',','precision',5);
dlmwrite('DA_spouseagediff_syndata1.txt',one_spouseagediff_syndata1','delimiter',',','precision',5);
dlmwrite('DA_spouseagediff_syndata2.txt',one_spouseagediff_syndata2','delimiter',',','precision',5);
dlmwrite('DA_spouseagediff_syndata3.txt',one_spouseagediff_syndata3','delimiter',',','precision',5);
dlmwrite('DA_spouseagediff_syndata4.txt',one_spouseagediff_syndata4','delimiter',',','precision',5);
dlmwrite('DA_spouseagediff_syndata5.txt',one_spouseagediff_syndata5','delimiter',',','precision',5);

dlmwrite('RS_spouseagediff_orig.txt',one_spouseagediff_orig','delimiter',',','precision',5);
dlmwrite('RS_spouseagediff_syndata1.txt',one_spouseagediff_syndata1','delimiter',',','precision',5);
dlmwrite('RS_spouseagediff_syndata2.txt',one_spouseagediff_syndata2','delimiter',',','precision',5);
dlmwrite('RS_spouseagediff_syndata3.txt',one_spouseagediff_syndata3','delimiter',',','precision',5);
dlmwrite('RS_spouseagediff_syndata4.txt',one_spouseagediff_syndata4','delimiter',',','precision',5);
dlmwrite('RS_spouseagediff_syndata5.txt',one_spouseagediff_syndata5','delimiter',',','precision',5);

%%
zHHunique = zeros(1,500);
zmemberunique = zeros(1,500);
for i=1:500
    zHHunique(i) = length(unique(results.z_HH_save((i-1)*10+5010,:)));
    zmemberunique(i) = length(unique(results.z_member_save((i-1)*10+5010,:)));
end

dlmwrite('ACSsimulation_nozeros_zHHunique_112514.txt',zHHunique,'delimiter',',','precision',5);
dlmwrite('ACSsimulation_nozeros_zmemberunique_112514.txt',zmemberunique,'delimiter',',','precision',5);
mean(zHHunique)
mean(zmemberunique)

mean(zHHunique)
quantile(zHHunique,.025)
quantile(zHHunique,.975)

min(zmemberunique)
max(zmemberunique)

nout = zeros(1,500);
for i=1:500
    nout(i) = results.nout_extend((i-1)*10+5010,:);
end
dlmwrite('ACSsimulation_nozeros_nout_042214.txt',nout,'delimiter',',','precision',7);

