%%%%---define output files---%%
output.alphaout = zeros(nrun,1); 
output.betaout = zeros(nrun,1);
output.piout = zeros(eff.sam,hyper.K); 
output.wout = zeros(eff.sam,hyper.K,hyper.L);
output.nout = zeros(nrun,1);
output.n_sout = zeros(nrun,1);

output.extrasize = zeros(nrun,3);

output.z_HH_save = zeros(nrun,orig.n_individuals);
output.z_member_save = zeros(nrun,orig.n_individuals);
output.elapsed_time = zeros(1,nrun);

output.newphiout = zeros(eff.sam,orig.maxd*orig.p,hyper.K*hyper.L);
output.lambda1out = zeros(eff.sam,hyper.K,hyper.dHH(1));
output.lambda2out = zeros(eff.sam,hyper.K,hyper.dHH(2));