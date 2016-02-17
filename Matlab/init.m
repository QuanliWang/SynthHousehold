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

initParameters;

initOutput;