addpath dist

n = 10000;
[SS,origdata] = PrepareData(n);

%%
% variables: HHindex 1, pernum 2, sex 3, race 4, ethn 5, age 6, relate 7(ind level) ownership 8 (binary, HH level)
[n_individuals,p] = size(origdata);
p = p-4;

HHrowIndex = [1 cumsum(SS)+1]; 
HHdataorig = origdata(HHrowIndex(1: (end -1)),8:9);
HHdata1_all = HHdataorig(:,1);
HHdata2_all = HHdataorig(:,2)-1;
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

%nrun=1000; burn=800; thin=10; %50;
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

extrasize = zeros(nrun,3);

z_HH_save = zeros(nrun,n_individuals);
z_member_save = zeros(nrun,n_individuals);
elapsed_time = zeros(1,nrun);

newphiout = zeros(eff.sam,maxd*p,K*L);
lambda1out = zeros(eff.sam,K,dHH1);
lambda2out = zeros(eff.sam,K,dHH2);

%%%%---initial values---%%
alpha=1;        % hyperparameters for stick-breaking weights
beta=1;

phi=zeros(maxd,p,K*L);      % cell probabilities
for j=1:p
    for c=1:d(j)
        phi(c,j,1)=sum(data(:,j)==c)/n_individuals;
    end
end
for count = 2:K*L
    phi(:,:,count) = phi(:,:,1);
end

%initialize lambda
lambda1 = zeros(K,dHH1);
for c=1:dHH1
    lambda1(:,c) = sum(HHdata1_all==c)/n;
end

lambda2 = zeros(K,dHH2);
for c=1:dHH2
    lambda2(:,c) = sum(HHdata2_all==c)/n;
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

%
data_full_all = origdata;
HHdata1_all = HHdataorig(:,1);
HHdata2_all = HHdataorig(:,2)-1;

ImpossibleIndividuals = [];z_HH_extra = [];