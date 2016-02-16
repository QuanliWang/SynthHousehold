addpath dist

%% initilize initial data
initInputData;

%%
%%%%---global parameters---%%
%nrun=10000; burn=8000; thin=10; %50;
%nrun=1000; burn=800; thin=10; %50;
nrun=25; burn=20; thin=1;
eff.sam=(nrun-burn)/thin;

initHyperParameters;

%%%%---initial values---%%
alpha=1;        % hyperparameters for stick-breaking weights
beta=1;

phi=zeros(orig.maxd,orig.p,hyper.K*hyper.L);      % cell probabilities
for j=1:orig.p
    for c=1:orig.d(j)
        phi(c,j,1)=sum(orig.data(:,j)==c)/orig.n_individuals;
    end
end
for count = 2:hyper.K*hyper.L
    phi(:,:,count) = phi(:,:,1);
end

HHdata_all = orig.HHdataorig; HHdata_all(:,2) = HHdata_all(:,2) - 1;

%initialize lambda
lambda1 = zeros(hyper.K,hyper.dHH(1));
for c=1:hyper.dHH(1)
    lambda1(:,c) = sum(HHdata_all(:,1)==c)/orig.n;
end

lambda2 = zeros(hyper.K,hyper.dHH(2));
for c=1:hyper.dHH(2)
    lambda2(:,c) = sum(HHdata_all(:,2)==c)/orig.n;
end

u=[betarnd(1,alpha,[hyper.K-1,1]);1];
pi=zeros(hyper.K,1);
pi(1:hyper.K)=u.*cumprod([1;1-u(1:hyper.K-1)]);

v=[betarnd(1,beta,[hyper.K,hyper.L-1]),ones(hyper.K,1)];
w=zeros(hyper.K,hyper.L);
v1 = (v(1,:))';
w(1,1:hyper.L)=v1.*cumprod([1;1-v1(1:hyper.L-1)]);
for k=1:hyper.K
    v1 = (v(k,:))';
    w(k,1:hyper.L)=v1.*cumprod([1;1-v1(1:hyper.L-1)]);
end

%
data_full_all = orig.origdata;

ImpossibleIndividuals = [];z_HH_extra = [];

initOutput;