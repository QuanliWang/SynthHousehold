addpath dist

%% initilize initial data
initInputData;

% what are those?
hyper.dHH = [2 3];

%%
%%%%---global parameters---%%
%nrun=10000; burn=8000; thin=10; %50;
%nrun=1000; burn=800; thin=10; %50;
nrun=25; burn=20; thin=1;
eff.sam=(nrun-burn)/thin;

K=40;                     
L=15;
aa=0.25; ab=0.25;                  % gamma hyperparameters for alpha
ba=0.25; bb=0.25;

%%%%---initial values---%%
alpha=1;        % hyperparameters for stick-breaking weights
beta=1;

phi=zeros(orig.maxd,orig.p,K*L);      % cell probabilities
for j=1:orig.p
    for c=1:orig.d(j)
        phi(c,j,1)=sum(orig.data(:,j)==c)/orig.n_individuals;
    end
end
for count = 2:K*L
    phi(:,:,count) = phi(:,:,1);
end

HHdata_all = orig.HHdataorig; HHdata_all(:,2) = HHdata_all(:,2) - 1;

%initialize lambda
lambda1 = zeros(K,dHH(1));
for c=1:dHH(1)
    lambda1(:,c) = sum(HHdata_all(:,1)==c)/orig.n;
end

lambda2 = zeros(K,dHH(2));
for c=1:dHH(2)
    lambda2(:,c) = sum(HHdata_all(:,2)==c)/orig.n;
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

%
data_full_all = orig.origdata;

ImpossibleIndividuals = [];z_HH_extra = [];


%%%%---define output files---%%
alphaout = zeros(nrun,1); 
betaout = zeros(nrun,1);
piout = zeros(eff.sam,K); 
wout = zeros(eff.sam,K,L);
nout = zeros(nrun,1);
n_sout = zeros(nrun,1);

extrasize = zeros(nrun,3);

z_HH_save = zeros(nrun,orig.n_individuals);
z_member_save = zeros(nrun,orig.n_individuals);
elapsed_time = zeros(1,nrun);

newphiout = zeros(eff.sam,orig.maxd*orig.p,K*L);
lambda1out = zeros(eff.sam,K,dHH(1));
lambda2out = zeros(eff.sam,K,dHH(2));