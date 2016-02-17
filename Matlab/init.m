addpath dist

%% initilize initial data
initInputData;

%%
%%%%---global parameters---%%
%nrun=10000; burn=8000; thin=10; %50;
%nrun=1000; burn=800; thin=10; %50;
mc.nrun=25; mc.burn=20; mc.thin=1;
mc.eff.sam=(mc.nrun-mc.burn)/mc.thin;

initHyperParameters;

%%%%---initial values---%%
para.alpha=1;        % hyperparameters for stick-breaking weights
para.beta=1;

para.phi=zeros(orig.maxd,orig.p,hyper.K*hyper.L);      % cell probabilities
for i=1:orig.p
    for j=1:orig.d(i)
        para.phi(j,i,1)=sum(orig.data(:,i)==j)/orig.n_individuals;
    end
end
for i = 2:hyper.K*hyper.L
    para.phi(:,:,i) = para.phi(:,:,1);
end

para.HHdata_all = orig.HHdataorig; para.HHdata_all(:,2) = para.HHdata_all(:,2) - 1;

%initialize lambda
para.lambda = cell(length(hyper.dHH),1);
for i = 1: length(hyper.dHH)
    lambda = zeros(hyper.K,hyper.dHH(i));
    for j=1:hyper.dHH(i)
        lambda(:,j) = sum(para.HHdata_all(:,i)==j)/orig.n;
    end
    para.lambda{i} = lambda;
end

para.u=[betarnd(1,para.alpha,[hyper.K-1,1]);1];
para.pi=zeros(hyper.K,1);
para.pi(1:hyper.K)=para.u.*cumprod([1;1-para.u(1:hyper.K-1)]);

para.v=[betarnd(1,para.beta,[hyper.K,hyper.L-1]),ones(hyper.K,1)];
para.w=zeros(hyper.K,hyper.L);

v1 = (para.v(1,:))';
para.w(1,1:hyper.L)=v1.*cumprod([1;1-v1(1:hyper.L-1)]);
for i=1:hyper.K
    v1 = (para.v(i,:))';
    para.w(i,1:hyper.L)=v1.*cumprod([1;1-v1(1:hyper.L-1)]);
end
clear v1

%
para.data_full_all = orig.origdata;

para.ImpossibleIndividuals = [];para.z_HH_extra = [];

initOutput;