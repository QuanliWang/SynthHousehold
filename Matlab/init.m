addpath dist

%% initilize initial data
initInputData;

%%
%%%%---global parameters---%%
%mc.nrun=10000; mc.burn=8000; mc.thin=10;   %long
%mc.nrun=1000; mc.burn=800; mc.thin=10;     %mediun
mc.nrun=25; mc.burn=20; mc.thin=1;          %short

mc.eff.sam=(mc.nrun-mc.burn)/mc.thin;

initHyperParameters;

initParameters;

initOutput;