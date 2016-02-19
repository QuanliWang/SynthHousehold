
if exist('household.mat', 'file') == 2
    disp('using existing data')
    load household
else
    disp('random sampling new data')
    n = 10000;
    [SS,origdata] = PrepareData('Monica/ACShouseholddata_size234all.mat',n);
    save household
    %dlmwrite('data.txt',origdata,'delimiter','\t');
end

%%
orig.origdata = origdata;
orig.SS = SS;
orig.n = length(orig.SS);
HHrowIndex = [1 cumsum(orig.SS)+1]; 
orig.HHdataorig = orig.origdata(HHrowIndex(1: (end -1)),8:9);
orig.HHserial = orig.origdata(:,1);
clear origdata SS n desc HHrowIndex

[orig.n_individuals,orig.p] = size(orig.origdata);
orig.p = orig.p-4;

% variables: HHindex 1, pernum 2, sex 3, race 4, ethn 5, age 6, relate
% 7(ind level) ownership 8 (binary, HH level) household size 
orig.d = [2,9,5,94,12]; %number of levels for each variable (5 variables all together)
orig.data = orig.origdata(:,3:7);
orig.dataT = orig.origdata(:,3:7)';
orig.maxd = max(orig.d);

orig.ACS_count = zeros(3,1);
orig.ACS_count(1) = sum(orig.SS==2);
orig.ACS_count(2) = sum(orig.SS==3);
orig.ACS_count(3) = sum(orig.SS==4);